#!/bin/bash -l

graph=$1
cur_dir=`pwd`

bench=SSSP
for bench in SSSP BFS SSWP SSNP Viterbi Heter; do
    mode=glign_intra
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    intra=`grep -e 'Glign-Intra evaluation time: ' $file | awk '{print $NF}' | bc`
    echo "${graph} ${bench} ${mode} Glign-Intra time: ${intra}s"

    mode=glign_inter
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    inter=`grep -e 'Glign-Inter evaluation time: ' $file | awk '{print $NF}'`
    speedup=$(echo "scale=2;(($intra / $inter))" | bc)
    echo "${graph} ${bench} ${mode} speedup: ${speedup}"
    echo "========"
done