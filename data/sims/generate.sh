#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <generate>" >&2
    exit 1
fi

for k in {3,4,5}
do
    for m in {2,3,4}
    do
	for n in {50,100,500}
	do
	    for seed in {1..20}
	    do
		echo Generating $k clusters, $m samples, $n SNVs, seed $seed
		$1 -C 200 -k $k -m $m -n $n -s $seed -maxCN 1 -o .
	    done
	done
    done
done
