import networkx as nx
import numpy as np
from depytrace.eval import conductance
from random import random
from collections import deque as Que


def eigen(L):
    u1 = np.ones((L.shape[1], 1))
    u1prev = 0
    while np.linalg.norm(u1 - u1prev, 1) > 1.E-3:
        u1prev = u1
        u1 = L @ u1
        u1 = u1 / np.linalg.norm(u1, ord=2)
    l1 = np.linalg.norm(L@u1, ord=2)/np.linalg.norm(u1, ord=2)
    return u1, l1


def laplacian(G, vertex2id):
    A = np.zeros((len(vertex2id), len(vertex2id)))
    for u, v in G.edges():
        if u != v:
            A[vertex2id[u]][vertex2id[v]] = 1
            A[vertex2id[v]][vertex2id[u]] = 1
    invD = np.diag((A.sum(axis=1) + 1.E-12) ** (-1))
    return np.eye(len(vertex2id)) - invD ** 0.5 @ A @ invD ** 0.5


def subgraph(G, sampler, r):
    T = nx.DiGraph([(u, v) for u, v in G.edges() if random() < sampler(u, v)])
    return bfs_tree(T, r)


def trivial_graph(r):
    T = nx.DiGraph()
    T.add_node(r)
    return T


def bfs_tree(T, r):
    if r not in T:
        return trivial_graph(r)
    return nx.traversal.bfs_tree(T, r)


def eigenreductor(G, r):
    vertex2id = {u: i for i, u in enumerate(G)}
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian(G, vertex2id))
    eigorder = sorted(list(range(len(eigenvalues))), key=lambda i: -eigenvalues[i])
    max_cond, max_tree = 0, None
    original_sum = eigenvalues.sum()
    largest = 0
    smallest = len(eigenvalues)-1
    for i in range(len(eigenvalues)-1):
        if eigenvalues[eigorder[smallest]] / (eigenvalues.sum() - eigenvalues[eigorder[largest]]) > eigenvalues[eigorder[smallest-1]]/(eigenvalues.sum() - eigenvalues[eigorder[smallest]]):
            changed = (eigorder[largest], eigenvalues[eigorder[largest]])
            eigenvalues[eigorder[largest]] = 0
            largest += 1
        else:
            changed = (eigorder[smallest], eigenvalues[eigorder[smallest]])
            eigenvalues[eigorder[smallest]] = 0
            smallest -= 1
        Lapprox = eigenvectors @ np.diag(eigenvalues*original_sum/eigenvalues.sum()) @ eigenvectors.transpose()
        #for Y in np.arange(0, 1, 0.1):
        subgraph = nx.DiGraph([(u,v) for u,v in G.edges() if random() < abs(Lapprox[vertex2id[u],vertex2id[v]])])
        if r not in subgraph:
            eigenvalues[changed[0]] = changed[1]
            continue
        tree = bfs_tree(subgraph, r)
        cond = conductance(G, tree)
        if cond > max_cond:
            max_cond, max_tree = cond, tree
    return trivial_graph(r) if max_tree is None else max_tree


class SOTASparsifier:
    def __init__(self):
        self.repeats = 10
        self.range_a = 10

    def __call__(self, G, r):
        # https://arxiv.org/pdf/0803.0929.pdf
        a = 0
        maxT = G
        maxConductance = 0
        vertex2id = {u: i for i, u in enumerate(G)}
        A = np.zeros((len(vertex2id), len(vertex2id)))
        for u, v in G.edges():
            if u != v:
                A[vertex2id[u]][vertex2id[v]] = 1
                A[vertex2id[v]][vertex2id[u]] = 1
        invD = np.diag((A.sum(axis=0) + 1.E-12) ** (-1))
        L = np.eye(len(vertex2id)) - invD ** 0.5 @ A @ invD ** 0.5
        pseudoinverse = 0
        eigenvalues, eigenvectors = np.linalg.eigh(L)
        for eigenvalue, eigenvector in zip(eigenvalues, eigenvectors):
            if eigenvalue > 1.E-12:
                pseudoinverse += np.outer(eigenvector,eigenvector)/eigenvalue

        epsilon = 0.1
        beta = 1  # 1+-epsilon factor approximation with probability 1-n**-beta (where n is the number of nodes)
        k = int((4+2*beta)*np.log(len(G))/(epsilon**2/2-epsilon**3/3))
        for repeat in range(10*self.range_a):
            if repeat % self.repeats == 0:
                a += 1
            Q = np.random.binomial(1, 0.5, (k, len(G)))*2-1
            Q = Q/(float(k)**0.5)
            Z = Q@pseudoinverse

            if random()<1./self.repeats:
                a += 1
            T = nx.DiGraph()
            Y = float(len(G))/sum(np.linalg.norm(Z[vertex2id[u],:]-Z[vertex2id[v],:])**2 for u, v in G.edges())/a
            for u, v in G.edges():
                resistance = np.linalg.norm(Z[vertex2id[u],:]-Z[vertex2id[v],:])**2
                if random() < Y/resistance:
                    T.add_edge(u, v)
            if r not in T:
                T = nx.DiGraph()
                T.add_node(r)
                continue
            else:
                T = nx.traversal.bfs_tree(T, r)
            cond = conductance(G, T)
            if cond > maxConductance:
                maxConductance = cond
                maxT = T
        #if maxT != G:
        #    return SOTASparsifier()(maxT, r)
        return maxT


class Greedy:
    def __init__(self, evaluator=lambda a,b: float(a-b)/b, threshold=1):
        self.evaluator = evaluator
        self.threshold = threshold

    def __call__(self, G, rs, Gconstrained=None):
        if Gconstrained is None:
            Gconstrained = G
        T = nx.DiGraph(rs.edges() if isinstance(rs, nx.DiGraph) else None)
        if not isinstance(rs, list):
            rs = [rs]
        for r in rs:
            T.add_node(r)
        visit = Que([r for r in rs])
        visited = set([r for r in rs])
        num_internal = rs.number_of_edges() if isinstance(rs, nx.DiGraph) else 0
        num_outgoing = sum(G.out_degree(u) for u in rs)
        while visit:
            u = visit.popleft()
            for v in Gconstrained.successors(u):
                if v in visited:
                    continue
                outgoing = [k for k in G.successors(v) if k not in T]
                if num_internal == 0 or self.evaluator(num_outgoing, num_internal) \
                        < self.threshold*self.evaluator(num_outgoing+len(outgoing), num_internal+1):
                    T.add_edge(u, v)
                    visit.append(v)
                    visited.add(v)
                    num_outgoing += len(outgoing)
                    num_internal += 1
        return T


def greedy_recursive(G, rs, evals):
    for eval in evals:
        rs = Greedy(eval)(G, rs)
    return rs


class Sparsifier:
    def __init__(self, measure=conductance, directed=True):
        self.measure = measure
        self.directed = directed

    def _directed_call(self, G, r):
        max_cond, max_tree = 0, None
        max_deg = max(G.out_degree(v) for v in G)
        for Y in np.arange(0, 1, 0.1):
            for _ in range(10):
                subgraph = nx.DiGraph()
                subgraph.add_node(r)
                next = Que([r])
                visited = set()
                while next:
                    v = next.popleft()
                    if v in visited:
                        continue
                    visited.add(v)
                    for u in G.successors(v):
                        if random() < Y*max_deg / G.out_degree(v):
                            subgraph.add_edge(v, u)
                            next.append(u)
                tree = bfs_tree(subgraph, r)
                cond = self.measure(G, tree)
                if r in subgraph and cond > max_cond:
                    max_cond, max_tree = cond, tree
        return trivial_graph(r) if max_tree is None else max_tree

    def _undirected_call(self, G, r):
            # https://arxiv.org/pdf/0808.4134.pdf
            max_cond, max_tree = 0, None
            max_deg = max(G.degree(v) for v in G)
            # print('max deg', max_deg)
            for Y in np.arange(0, 1, 0.01):
                for _ in range(1):
                    subgraph = nx.DiGraph(
                        [(u, v) for u, v in G.edges() if random() < Y * max_deg / min(G.degree(u), G.degree(v))])
                    # tree = subgraph
                    tree = bfs_tree(subgraph, r)
                    cond = self.measure(G, tree)
                    if cond > max_cond:
                        max_cond, max_tree = cond, tree
            return trivial_graph(r) if max_tree is None else max_tree

    def __call__(self, G, r):
        # https://arxiv.org/pdf/0808.4134.pdf
        return self._directed_call(G, r) if self.directed else self._undirected_call(G, r)
