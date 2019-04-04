#!/bin/bash
if [ ! $# -eq 1 ]
then
    echo "Usage: $0 <generate>" >&2
    exit 1
fi

for k in {8,10,15}
do
    for m in {4,6,8}
    do
	for n in {100,500}
	do
	    for seed in {1..20}
	    do
		echo Generating $k clusters, $m samples, $n SNVs, seed $seed
		$1 -C 200 -k $k -m $m -n $n -s $seed -maxCN 1 -o .
	    done
	done
    done
done
