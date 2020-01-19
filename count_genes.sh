#!/bin/bash
input=$1


total_counts=`awk '{ print }' $input | wc -l`
uniq_gene_counts=`awk '!x[$0]++' $input | wc -l`
dup_gene_couns=`awk 'x[$0]++' $input | awk '!x[$0]++' | wc -l`

cat << EOF
Total Counts: $total_counts
Unique Counts: $uniq_gene_counts
Duplicated Counts: $dup_gene_couns
EOF