#!/bin/bash



echo "method	samples	SNVs	clusters	coverage	seed	ARI	RI	dist	recall	precision	inferred_clusters"
for m in {2,3,4}
do
    for k in {3,4,5}
    do
	for F in ../results/sims/m${m}*k${k}*
	do
	    echo -n "DECIFER	"
	    BASE=$(basename ${F})

	    BEST="$(python - $F <<END
import sys, os, glob

b = os.path.basename(sys.argv[1])
with open(os.path.join(sys.argv[1], '{}.{}'.format(b, 'BIC.tsv')), 'r') as i:
    i.readline()
    bic = {int(l.strip().split()[0]) : float(l.strip().split()[-1]) for l in i if len(l) > 1}

print min(bic.keys(), key=(lambda x : bic[x]))

END
)"

	    python processBIC.py $k $m ../data/sims/${BASE}.SNV.tsv ${F}/${BASE}_k${BEST}.SNV.tsv
	done
    done
done
