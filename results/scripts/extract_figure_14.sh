#!/bin/bash -l

graph=$1
cur_dir=`pwd`

bench=SSSP
for bench in SSSP BFS SSWP SSNP Viterbi Heter; do
    total_F=0;
    mode=ligra_s
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    total_F=`grep -e 'Sequential F: ' $file | awk '{SUM += $NF} END { print SUM }'`
    # echo "${graph} ${bench} ${mode} seq F: ${total_F}"

    mode=glign_intra
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    intra_F=`grep -e 'Glign-Intra F: ' $file | awk '{SUM += $NF} END { print SUM }'`
    # echo "${graph} ${bench} ${mode} intra F: ${intra_F}"
    affinity=$(echo "scale=4;(($intra_F / $total_F))" | bc)
    echo "${graph} ${bench} ${mode} 1-Affinity: ${affinity}"

    mode=glign_inter
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    inter_F=`grep -e 'Glign-Inter F: ' $file | awk '{SUM += $NF} END { print SUM }'`
    # echo "${graph} ${bench} ${mode} inter F: ${inter_F}"
    affinity=$(echo "scale=4;(($inter_F / $total_F))" | bc)
    echo "${graph} ${bench} ${mode} 1-Affinity: ${affinity}"

    mode=glign_batch
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    batch_F=`grep -e 'Glign-Batch F: ' $file | awk '{SUM += $NF} END { print SUM }'`
    # echo "${graph} ${bench} ${mode} intra F: ${batch_F}"
    affinity=$(echo "scale=4;(($batch_F / $total_F))" | bc)
    echo "${graph} ${bench} ${mode} 1-Affinity: ${affinity}"

    # mode=glign
    # file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    # glign=`grep -e 'Glign evaluation time: ' $file | awk '{print $NF}'`
    # speedup=$(echo "scale=2;(($seq / $glign))" | bc)
    # echo "${graph} ${bench} ${mode} speedup: ${speedup}"
    echo "========"
done