#!/bin/bash

# for each */genes-regulate-genes.txt, get the list of gene names

for f in *
do
    if [[ -d "$f" ]]
    then
	fn="${f}/genes-regulate-genes.txt"
	echo "${fn}"
	cat <(cut -f1 "$fn") <(cut -f2 "$fn") | sort -u > "${f}/genes.txt"
    fi
done

