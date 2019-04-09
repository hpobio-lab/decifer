#!/usr/bin/python
import clustering as clust
import sys
import networkx as nx

def parse_true_sol(k, m, filename):
    C = [ [] for j in range(k) ]
    dcf = [ [ 0 for j in range(k) ] for p in range(m) ]

    with open(filename) as f:
        f.readline()
        for line in f:
            s = line.rstrip("\n").split("\t")
            SNV = int(s[0])
            cluster = int(s[1])
            C[cluster].append(SNV)
            for p in range(m):
                dcf[p][cluster] = float(s[4 + 3 * p])

    return C, dcf

def parse_inf_sol(k, m, filename):
    C = [ [] for j in range(k) ]
    dcf = [ [ 0 for j in range(k) ] for p in range(m) ]

    with open(filename) as f:
        f.readline()
        for line in f:
            s = line.rstrip("\n").split("\t")
            SNV = int(s[0])
            cluster = int(s[k+1])
            C[cluster].append(SNV)
            for p in range(m):
                dcf[p][cluster] = float(s[k+2+p])

    return C, dcf

def readk(inf_sol_filename):
    with open(inf_sol_filename, 'r') as i:
        hd = i.readline().strip().split()
        return sum('cluster' in e for e in hd) - 1

def compare_dcfs(k, m, dcf1, dcf2):
    assert len(dcf1[0]) == k
    k2 = len(dcf2[0])

    tot_dist = 0
    for j1 in range(k):
        min_d = 1e100
        for j2 in range(k2):
            d = 0.
            for p in range(m):
                d += abs(dcf1[p][j1] - dcf2[p][j2-k2])
            if d < min_d:
                min_d = d
        tot_dist += min_d

    for j2 in range(k2):
        min_d = 1e100
        for j1 in range(k):
            d = 0.
            for p in range(m):
                d += abs(dcf1[p][j1] - dcf2[p][j2-k2])
            if d < min_d:
                min_d = d
        tot_dist += min_d

    return tot_dist / ((k + k2) * m)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.stderr.write("Usage: %s <k> <m> <true_solution> <inferred_solution>\n" % sys.argv[0])
        sys.exit(1)

    k_true = int(sys.argv[1])
    m = int(sys.argv[2])
    true_sol_filename = sys.argv[3].replace('_BEST', '')
    inf_sol_filename = sys.argv[4]

    C_true, DCF_true = parse_true_sol(k_true, m, true_sol_filename)

    k_inf = readk(inf_sol_filename)
    C_inf, DCF_inf = parse_inf_sol(k_inf, m, inf_sol_filename)

    ari = clust.adjusted_rand_index(C_true, C_inf)
    ri = clust.rand_index(C_true, C_inf)
    dist = compare_dcfs(k_true, m, DCF_true, DCF_inf)
    recall, precision = clust.recall_and_precision(C_inf, C_true)
    s = true_sol_filename.split("/")[-1].rstrip(".SNV.tsv").split("_")

    print "\t".join([s[0][1:], s[1][1:], s[2][1:], s[3][1:], s[4][1:], str(ari), str(ri), str(dist), str(recall), str(precision), str(k_inf)])
