#!/bin/bash -l

graph=$1
cur_dir=`pwd`

bench=SSSP
for bench in SSSP BFS SSWP SSNP Viterbi Heter; do
    mode=glign_intra
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    intra=`grep -e 'Glign-Intra evaluation time: ' $file | awk '{print $NF}' | bc`
    echo "${graph} ${bench} ${mode} Glign-Intra time: ${intra}s"

    mode=glign_batch
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    glign=`grep -e 'Glign-Batch evaluation time: ' $file | awk '{print $NF}'`
    speedup=$(echo "scale=2;(($intra / $glign))" | bc)
    echo "${graph} ${bench} ${mode} speedup: ${speedup}"
    echo "========"
done