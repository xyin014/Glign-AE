#!/bin/bash -l

graph=$1
cur_dir=`pwd`

bench=BFS

mode=ibfs
file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
bfs=`grep -e 'Sequential evaluation time: ' $file | awk '{print $NF}' | bc`
echo "${graph} ${bench} ${mode} iBFS time: ${bfs}s"

mode=glign_intra
file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
intra=`grep -e 'Glign-Intra evaluation time: ' $file | awk '{print $NF}'`
speedup=$(echo "scale=2;(($bfs / $intra))" | bc)
echo "${graph} ${bench} ${mode} speedup: ${speedup}"

mode=glign_batch
file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
batch=`grep -e 'Glign-Batch evaluation time: ' $file | awk '{print $NF}'`
speedup=$(echo "scale=2;(($bfs / $batch))" | bc)
echo "${graph} ${bench} ${mode} speedup: ${speedup}"