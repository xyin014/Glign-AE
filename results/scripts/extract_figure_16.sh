#!/bin/bash -l

graph=$1
cur_dir=`pwd`

bench=SSSP
for bench in SSSP BFS SSWP Viterbi; do
    
    mode=ligra_s
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    seq=`grep -e 'Sequential evaluation time: ' $file | awk '{print $NF}' | bc`
    echo "${size} ${graph} ${bench} ${mode} seq time: ${seq}s"

    for size in 2 4 8 16 32 64 128; do
        mode=ligra_c
        file="${cur_dir}/${graph}/${bench}_${mode}_${graph}_${size}.txt"
        ligra_c=`grep -e 'Ligra-C time: ' $file | awk '{print $NF}' | bc`
        speedup=$(echo "scale=2;(($seq / $ligra_c))" | bc)
        echo "Batch size ${size} ${graph} ${bench} ${mode} speedup: ${speedup}"

        mode=glign_batch
        file="${cur_dir}/${graph}/${bench}_${mode}_${graph}_${size}.txt"
        glign=`grep -e 'Glign-Batch evaluation time: ' $file | awk '{print $NF}'`
        speedup=$(echo "scale=2;(($seq / $glign))" | bc)
        echo "Batch size ${size} ${graph} ${bench} ${mode} speedup: ${speedup}"

        mode=glign_inter
        file="${cur_dir}/${graph}/${bench}_${mode}_${graph}_${size}.txt"
        glign=`grep -e 'Glign-Inter evaluation time: ' $file | awk '{print $NF}'`
        speedup=$(echo "scale=2;(($seq / $glign))" | bc)
        echo "Batch size ${size} ${graph} ${bench} ${mode} speedup: ${speedup}"
    done
    echo "========"
done