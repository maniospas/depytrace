import trace
import eval
from random import random, seed
import glob
from timeit import default_timer as time
from matplotlib import pyplot as plt


def er_graphs():
    for i in range(40):
        seed(i)
        nodes = int(100+random()*400)
        graph_type = ('ER', nodes, 0.01+random()*10/nodes)
        yield eval.create_graph(*graph_type, i, one_direction=False, random_weights=int(random()*10))


def dependeny_graphs():
    datasets = {file.replace("_", "")[len("data/dependencies") + 1:-4]: file for file in
                glob.glob('data/dependencies/*.csv') if "package" not in file}
    for dataset in datasets:
        yield eval.import_graph(datasets[dataset], directed=False).to_directed()


implementation = {
    #"eigen": trace.eigenreductor,
    "greedy": trace.Greedy(),
    "sparse": trace.Sparsifier(eval.ELOD),
    "maxcut": trace.RPCST,
    "rpcst": trace.Core(trace.RPCST),
    "rpcst*": trace.Core(trace.cleverRPCST),
    #"fast": trace.Core(trace.ELODFast())
}

nodes = [len(G) for G in dependeny_graphs()]
edges = [G.number_of_edges() for G in dependeny_graphs()]
plt.subplot(3, 1, 1)
plt.bar(range(len(nodes)), [node for _, node in sorted(zip(edges, nodes))])
plt.yscale('log')
plt.ylabel('Nodes')
plt.gca().axes.xaxis.set_ticklabels([])
plt.subplot(3, 1, 2)
plt.bar(range(len(edges)), sorted(edges))
plt.yscale('log')
plt.ylabel('Edges')
plt.gca().axes.xaxis.set_ticklabels([])
plt.subplot(3, 1, 3)
plt.bar(range(len(nodes)), [edge/node for edge, node in sorted(zip(edges, nodes))])
plt.ylabel('Degrees')
plt.gca().axes.xaxis.set_ticklabels([])
plt.xlabel('Dataset')
plt.show()

conductances = {method: list() for method in implementation}
heterogenity = {method: list() for method in implementation}
edges = list()
times = {method: list() for method in implementation}

for G in dependeny_graphs():
    nodes = [v for v in G if G.out_degree[v] >= 2]
    if nodes:
        for _ in range(10):
            r = nodes[(int)(random()*len(nodes))]
            edges.append(G.number_of_edges())
            for method in implementation:
                tic = time()
                subgraph = implementation[method](G, r)
                times[method].append(time()-tic)
                conductances[method].append(eval.conductance(G, subgraph))
                heterogenity[method].append(eval.heterogenity(G, subgraph))

    R, crit = eval.friedman_ranks(conductances)
    print(f'Critical difference {crit:.1f}')
    for method in implementation:
        print(f"{method} \t {sum(conductances[method])/len(conductances[method]):.0f} ({R[method]:.1f}) in {sum(times[method])/len(times[method])*1000:.0f} ms,"
              f"{sum(heterogenity[method])/len(heterogenity[method]):.2f}")

from numpy.polynomial.polynomial import polyfit
import numpy as np
for i, method in enumerate(times):
    plt.subplot(5, 1, i+1)
    b, m = polyfit(edges, np.array(times[method])*1000, 1)
    plt.scatter(edges, np.array(times[method])*1000, label=method)
    plt.plot(sorted(edges), np.array(sorted(edges))*m+b, 'r-')
    #plt.legend()
    if i != len(times)-1:
        plt.gca().axes.xaxis.set_ticklabels([])
    plt.title(method)
    plt.ylabel("Time (ms)")
plt.xlabel("Edges")
plt.show()

