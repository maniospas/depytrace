
def spread2(G, subgraph, r):
    speed = 0
    p = {v: 1 if v == r else 0 for v in G}
    d = {v: 2 for v in G}
    for v in subgraph:
        d[v] = 0.1
    p0 = p
    for speed in range(100):
        p_next = {v: sum(p[u]/G.out_degree(u)*p[u]*d[u] for u in G._pred[v])*0.9+0.1*p0[v] for v in G}
        p_next_sum = sum(p.values())
        p_next = {v: p_next[v]/p_next_sum for v in G}
        if sum(abs(p[v]-p_next[v]) for v in G) < 1.E-6*len(G):
            break
        p = p_next
    return speed

from math import log

def spread(G, subgraph, r):
    p = {v: 0 for v in G}
    susceptible = {v: 1 for v in G}
    removed = {v: 0 for v in G}
    p[r] = 1
    susceptible[r] = 0

    d = {v: 1 for v in G}
    for v in subgraph:
        d[v] = float(len(subgraph))/len(G)
    #d[r] = 0.1

    p[r] = 1
    beta = 0.350
    gamma = 0.035 #recovery rate
    dt = 0.4#04
    max_infected = 0
    for speed in range(200):
        p_next = {v: p[v]+dt*(sum(p[u]*d[u]/G.out_degree(u) for u in G._pred[v])*susceptible[v]*beta-gamma*p[v]) for v in G}
        susceptible = {v: susceptible[v]-dt*beta*susceptible[v]*sum(p[u]*d[u]/G.out_degree(u) for u in G._pred[v]) for v in G}
        removed = {v: removed[v]+dt*gamma*p[v] for v in G}
        max_infected = max(max_infected, sum(p_next.values()))
        p = p_next
    #print(max_infected, len(subgraph))
    #infected = sum(p_next.values())/total*len(p_next)
    return -max_infected
