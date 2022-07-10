import networkx as nx


def _count(G, subgraph):
    internal = 0
    outgoing_total = 0
    for v in subgraph:
        internal += subgraph.out_degree[v]
        outgoing_total += G.out_degree[v]
    return internal, outgoing_total


def conductance(G, subgraph):
    if subgraph is None:
        raise Exception("None subgraph not valid")
    if not nx.is_connected(subgraph.to_undirected()):
        raise Exception("Non-connected subgraph")
    internal = 0
    outgoing_total = 0
    for v in subgraph:
        internal += subgraph.out_degree[v]
        outgoing_total += G.out_degree[v]
    if internal == 0:
        return 0
    return outgoing_total/internal-1


def ELOD(G, subgraph, a=1):
    internal = 0
    outgoing_total = 0
    for v in subgraph:
        internal += subgraph.out_degree[v]
        outgoing_total += G.out_degree[v]
    return outgoing_total-a*internal

#evaluate GW https://arxiv.org/pdf/1908.05767.pdf, optimal value 0.7632


def heterogenity(G, subgraph):
    if len(subgraph) <= 1:
        return 0
    packages = set([".".join(v.split(".")[:-1]) for v in subgraph])
    return (len(packages)-1)/(len(subgraph)-1)
