#!/bin/bash -l

graph=$1
cur_dir=`pwd`

bench=SSSP
for bench in SSSP BFS SSWP SSNP Viterbi Heter; do

    mode=ligra_c
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    if [ "$bench" == "Heter" ]; then 
        ligra_c=`grep -e 'Ligra-C evaluation time: ' $file | awk '{print $NF}' | bc`
    else 
        ligrac=`grep -e 'Ligra-C time: ' $file | awk '{print $NF}' | bc`
    fi
    echo "${graph} ${bench} ${mode} Ligra-C time: ${ligrac}s"

    mode=glign_intra
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    intra=`grep -e 'Glign-Intra evaluation time: ' $file | awk '{print $NF}'`
    speedup=$(echo "scale=2;(($ligrac / $intra))" | bc)
    echo "${graph} ${bench} ${mode} speedup: ${speedup}"
    echo "========"
done