import numpy as np
from collections import deque as Que
import networkx as nx
from .backup import ELODFast


def ElODTree(G, r, a, D0):
    def _backward(self, T, a, D0):
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

    def _forward(self, T, r, a, D):
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

    if D0 is None:
        D0 = {v: G.out_degree[v] for v in G}
    D = _backward(G, a, D0)
    return _forward(G, r, a, D)


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


def cleverRPCST(G, r, a, D=None, pruning='strong'):
    import pcst_fast
    node_map = {}
    node_map_inv = {}
    for v in G:
        node_map_inv[len(node_map)] = v
        node_map[v] = len(node_map)
    if D is None:
        D = {v: G.out_degree(v) for v in G}
    #D[r] = sum(D[v] for v in G)
    edges = list()
    for i, j in G.edges():
        edges.append([node_map[i], node_map[j]])
    if a < 0:
        a = 0
    vertices_selected, edges_selected = pcst_fast.pcst_fast(np.asarray(edges, dtype=np.int64),
                                   np.asarray([D[v]-min(a,D[v]) for v in G], np.float64),
                                   np.asarray([a-min(a,D[v]) for u, v in G.edges()], np.float64),
                                   node_map[r], 1, pruning, 0)
    T = nx.DiGraph()
    #print([(node_map_inv[edges[edge][0]], node_map_inv[edges[edge][1]]) for edge in edges_selected])
    for edge in edges_selected:
        if G.has_edge(node_map_inv[edges[edge][0]], node_map_inv[edges[edge][1]]):
            T.add_edge(node_map_inv[edges[edge][0]], node_map_inv[edges[edge][1]])
        if G.has_edge(node_map_inv[edges[edge][1]], node_map_inv[edges[edge][0]]):
            T.add_edge(node_map_inv[edges[edge][1]], node_map_inv[edges[edge][0]])
    if len(T) == 0:
        T = nx.DiGraph()
        T.add_node(r)
        return T
    T = nx.traversal.bfs_tree(T, r)
    return T


def RPCST(G, r, a=0, D=None, pruning='strong'):
    import pcst_fast
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
                                   node_map[r], 1, pruning, 0)
    T = nx.DiGraph()
    for edge in edges_selected:
        T.add_edge(node_map_inv[edges[edge][0]], node_map_inv[edges[edge][1]])
        #T.add_edge(node_map_inv[edges[edge][1]], node_map_inv[edges[edge][0]])
    if len(T) == 0:
        T = nx.DiGraph()
        T.add_node(r)
        return T
    T = nx.traversal.bfs_tree(T, r)
    return T


def maxELOD (G, r, a):
    pass


class Core:
    def __init__(self, trace_method=ELODFast(), eps=1.E-6, starting_a=-1):
        self.trace_method = trace_method
        self.eps = eps
        self.starting_a = starting_a

    def __call__(self, G, r):
        min_a = self.starting_a
        max_a = max(G.degree(v) for v in G)
        while abs(min_a - max_a) > self.eps:
            mid = (max_a + min_a) * 0.5
            trace = self.trace_method(G, r, mid)
            if trace is None or len(trace) <= 1:
                max_a = mid
            else:
                min_a = mid
        return self.trace_method(G, r, min_a)
