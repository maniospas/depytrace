import networkx as nx
import pcst_fast
import numpy as np
from depytrace.trace.utils.heap import Heap
from collections import deque as Que

def _backward(T, a, D0):
    D = {v: D0.get(v, 0) for v in T}
    need_calc = {v: T.out_degree(v) for v in T}
    pending = [v for v in T if need_calc[v] == 0]
    while pending:
        v = pending.pop(0)
        if T.in_degree(v) == 0:
            continue
        if T.in_degree(v) > 1:
            raise Exception("Can only perform the BACKWARD step in trees")
        p = next(T.predecessors(v))
        D[p] += max(D[v]-a, 0)
        need_calc[p] -= 1
        if need_calc[p] == 0:
            pending.append(p)
    return D


def _forward(T, r, a, D):
    trace = nx.DiGraph()
    trace.add_node(r)
    pending = list(T.successors(r))
    while pending:
        v = pending.pop(0)
        if T.in_degree(v) > 1:
            raise Exception("Can only perform the FORWARD step in trees")
        if D[v] >= a:
            trace.add_edge(next(T.predecessors(v)), v)
            pending.extend(T.successors(v))
    return trace


def _remove_cycles(G, r):
    a = max(deg for _, deg in G.out_degree())
    subtraceELOD = {v: float('-inf') for v in G}
    subtraceELOD[r] = G.out_degree[r]
    heap = Heap()
    visited = {v: False for v in G}
    subtraceELOD[r] = G.out_degree[r]
    heap.add(r, subtraceELOD[r])
    T = nx.DiGraph()
    T.add_node(r)
    for u in heap:
        visited[u] = True
        for v in G.successors(u):
            if not visited[v]:
                subtraceELOD[v] = max(subtraceELOD[v], G.out_degree[v]+subtraceELOD[u]-a)
                heap.add(v, -subtraceELOD[v])
                T.add_edge(u, v)
    return T


def _add(dict1, dict2):
    for u, value in dict2.items():
        dict1[u] = dict1.get(u, 0) + dict2[u]


def _traverse_order(G, r):
    in_degree = {v: G.in_degree[v] for v in G}
    should_visit = {v: in_degree[v] for v in G}
    pending = Que([r])
    pending_multi_in = Que()
    visit_order = list()
    while True:
        if pending:
            u = pending.popleft()
        elif pending_multi_in:
            u = pending_multi_in.popleft()
        else:
            break
        visit_order.append(u)
        for v in G.successors(u):
            should_visit[v] -= 1
            if should_visit[v] == 0:
                if in_degree[v] > 1:
                    pending_multi_in.append(v)
                else:
                    pending.append(v)
    if len(visit_order) != len(G):
        raise Exception("Topological order is defined only for connected graphs")
    return visit_order


class ELODTree:
    def __call__(self, G, r, a):
        D0 = {v: G.out_degree(v) for v in G}
        D = _backward(G, a, D0)
        return _forward(G, r, a, D)


class ELODFast:
    def __init__(self, strong=False):
        self.strong = strong
    """def _forward_dag(self, G, r, a, D):
        trace = nx.DiGraph()
        trace.add_node(r)
        pending = Heap()
        pending.add(r, 0)
        visited_from = {v: None for v in D}
        for u in pending:
            for v in G.successors(u):
                Dv = D[u][v]
                if Dv >= a and visited_from[v] is None:
                    pending.add(v, -Dv)
                    trace.add_edge(u, v)
                    visited_from[v] = u
        return trace
    """

    def _forward_dag(self, G, r, a, D):
        trace = nx.DiGraph()
        trace.add_node(r)
        pending = Heap()
        pending.add(r, 0)
        visited_from = {v: None for v in D}
        for u in pending:
            for v in G.successors(u):
                Dv = D[v]#Duv[u].get(v,0)
                if Dv >= a and visited_from[v] is None:
                    pending.add(v, Dv)
                    trace.add_edge(u, v)
                    visited_from[v] = u
        return trace

    def _backward_dag(self, G, r, a, D0, from_id=0, to_id=None, strong=False):
        accounted_successors = {v: set() for v in G}
        if to_id is None:
            to_id = len(self.order)
        order = self.order[from_id:to_id]
        order.reverse()
        D = {v: D0.get(v, 0) for v in G}
        Duv = {v: dict() for v in G}
        for p in order:
            if D[p] == 0:
                continue
            for v in G.successors(p):#sorted(G.successors(p), key=lambda v: -D[v]):
                Dv = D[v]
                if strong:
                    common_accounted = set(s for s in accounted_successors[p] if s in accounted_successors[v])
                    for s in sorted(common_accounted):
                        if s not in common_accounted:
                            continue
                        Dv -= max(D[s] - a, 0)
                        common_accounted -= accounted_successors[s]

                # check if the node should be added
                if Dv >= a:
                    Duv[p][v] = Dv
                    if strong:
                        for s in accounted_successors[v]:
                            accounted_successors[p].add(s)
                        if G.in_degree[v] > 1:
                            accounted_successors[p].add(v)
                    D[p] += Dv - a
        return D

    def __call__(self, G, r, a):
        D0 = {v: G.out_degree[v] for v in G}
        G = _remove_cycles(G, r)

        self.order = _traverse_order(G, r)
        self.node_id = {v: i for i, v in enumerate(self.order)}

        D = self._backward_dag(G, r, a, D0, self.strong)
        trace = self._forward_dag(G, r, a, D)
        if self.strong:
            D = _backward(trace, a, D0)
            trace = _forward(trace, r, a, D)
            return trace
        if len(trace)==1:
            T = nx.DiGraph()
            T.add_node(r)
            return T
        return PCSTFast()(trace, r, a, D0)


class MinCut:
    def __call__(self, G, r):
        from pygrank.algorithms.pagerank import PageRank
        G = G.to_directed()
        ranks = PageRank().rank(G, {r: 1})
        ranks = {v: ranks[v]/G.degree(v) for v in G}
        max_grap = 0
        threshold = 0
        prev_rank = None
        for v, rank in sorted(ranks.items(), key=lambda item: item[1], reverse=True):
            if prev_rank is not None:
                gap = (prev_rank-rank)
                print(gap)
                if gap > max_grap:
                    max_grap = gap
                    threshold = rank
            prev_rank = rank
        T = nx.DiGraph()
        T.add_node(r)
        for u, v in G.edges():
            if ranks[u] <= threshold and ranks[v] <= threshold:
                T.add_edge(u, v)
        return nx.ego_graph(T, r, radius=1000000)


class PCSTFast:
    def __call__(self, G, r, a, D=None):
        node_map = {}
        node_map_inv = {}
        for v in G:
            node_map_inv[len(node_map)] = v
            node_map[v] = len(node_map)
        #max_degree = max(G.out_degree(v) for v in G)
        if D is None:
            D = {v: G.out_degree(v) for v in G}
        edges = list()
        for i, j in G.edges():
            edges.append([node_map[i], node_map[j]])
        if a < 0:
            a = 0
        vertices_selected, edges_selected = pcst_fast.pcst_fast(np.asarray(edges, dtype=np.int64),
                                       np.asarray([D[v] for v in G], np.float64),
                                       np.asarray([a for _ in range(len(edges))], np.float64),
                                       node_map[r], 1, 'strong', 0)
        T = nx.Graph()
        for edge in edges_selected:
            T.add_edge(node_map_inv[edges[edge][0]], node_map_inv[edges[edge][1]])
        if len(T) == 0:
            T = nx.DiGraph()
            T.add_node(r)
            return T
        T = nx.traversal.bfs_tree(T, r)
        return T