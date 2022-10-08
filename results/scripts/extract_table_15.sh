#!/bin/bash -l

graph=$1
cur_dir=`pwd`

bench=SSSP
for bench in SSSP BFS SSWP; do
    current_bench=$bench
    mode=ligra_s
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    seq=`grep -e 'Sequential evaluation time: ' $file | awk '{print $NF}' | bc`
    echo "${graph} ${bench} ${mode} seq time: ${seq}s"

    mode=ligra_c
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    if [ "$bench" == "Heter" ]; then 
        ligra_c=`grep -e 'Ligra-C evaluation time: ' $file | awk '{print $NF}' | bc`
    else 
        ligra_c=`grep -e 'Ligra-C time: ' $file | awk '{print $NF}' | bc`
    fi
    speedup=$(echo "scale=2;(($seq / $ligra_c))" | bc)
    echo "${graph} ${bench} ${mode} speedup: ${speedup}"

    mode=glign_intra
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    intra=`grep -e 'Glign-Intra evaluation time: ' $file | awk '{print $NF}'`
    speedup=$(echo "scale=2;(($seq / $intra))" | bc)
    echo "${graph} ${bench} ${mode} speedup: ${speedup}"

    mode=glign_inter
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    inter=`grep -e 'Glign-Inter evaluation time: ' $file | awk '{print $NF}'`
    speedup=$(echo "scale=2;(($seq / $inter))" | bc)
    echo "${graph} ${bench} ${mode} speedup: ${speedup}"

    mode=glign_batch
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    glign=`grep -e 'Glign-Batch evaluation time: ' $file | awk '{print $NF}'`
    speedup=$(echo "scale=2;(($seq / $glign))" | bc)
    echo "${graph} ${bench} ${mode} speedup: ${speedup}"

    mode=glign
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    glign=`grep -e 'Glign evaluation time: ' $file | awk '{print $NF}'`
    speedup=$(echo "scale=2;(($seq / $glign))" | bc)
    echo "${graph} ${bench} ${mode} speedup: ${speedup}"

    echo "========"
done