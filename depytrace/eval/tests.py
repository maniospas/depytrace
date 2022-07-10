import scipy.stats

# https://www.real-statistics.com/statistics-tables/studentized-range-q-table/ for df=inf
studentizedq = {0.01: [None, 3.643, 4.120, 3.403, 4.603, 4.757, 4.882, 4.987, 5.078, 5.157, 5.227, 5.290]}


def compare(list1, list2):
    return str(sum([1. for u, v in zip(list1, list2) if u>v])/len(list1))+" ("+str(sum([1. for u, v in zip(list1, list2) if u==v])/len(list1))+" ties)"


def friedman_ranks(results):
    n = len(results[list(results.keys())[0]])
    total_ranks = {method: list() for method in results}
    for i in range(n):
        for method, rank in zip(results.keys(), scipy.stats.rankdata([-results[method][i] for method in results])):
            total_ranks[method].append(rank)
    se = (len(results)*(len(results)+1)/n/12)**0.5
    qcrit = studentizedq[0.01][len(results)]
    R = {method: sum(total_ranks[method])/n for method in results}
    return R, qcrit*se
