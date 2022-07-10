import networkx as nx
from random import random


def create_graph(method, nodes, param, seed, one_direction=False, random_weights=15):
    if method == 'ER':
        G = nx.erdos_renyi_graph(nodes, param, seed=seed).to_directed()
    elif method == 'BA':
        G = nx.generators.barabasi_albert_graph(nodes, param, seed=seed).to_directed()
    else:
        raise Exception("Invalid method")
    G.remove_nodes_from(list(nx.isolates(G)))
    G.remove_edges_from(nx.selfloop_edges(G))
    #G = max_conductance._remove_cycles(G, r)
    if one_direction:
        for u, v in list(G.edges()):
            if G.has_edge(u,v) and G.has_edge(v, u):
                G.remove_edge(v, u)
    extra_node = 0
    for u in list(G.nodes()):
        for _ in range((int)(random()*random_weights)):
            G.add_edge(u, "extra"+str(extra_node))
            if not one_direction:
                G.add_edge("extra"+str(extra_node), u)
            extra_node += 1
    return G


def import_graph(path, directed=False, delim=";", edge_cols=(0, 1), skip_title=False, postprocess=lambda x: x):
    G = nx.DiGraph() if directed else nx.Graph()
    with open(path) as file:
        if skip_title:
            next(file)
        for line in file:
            splt = line[:-1].split(delim)
            if len(splt) >= 2:
                G.add_edge(postprocess(splt[edge_cols[0]]), postprocess(splt[edge_cols[1]]))
    return G
