#!/bin/bash -l

graph=$1
cur_dir=`pwd`

bench=SSSP
for bench in SSSP BFS SSWP SSNP Viterbi Heter; do
    
    mode=glign_intra
    file="${cur_dir}/perf/perf_${bench}_${mode}_${graph}.log"
    tmp=`grep -e 'LLC-load-misses:u' $file`
    intra_cnt=`echo ${tmp%'LLC-load-misses:u'*}`

    mode=glign_inter
    file="${cur_dir}/perf/perf_${bench}_${mode}_${graph}.log"
    tmp=`grep -e 'LLC-load-misses:u' $file`
    inter_cnt=`echo ${tmp%'LLC-load-misses:u'*}`
    reduction=$(echo "scale=2;((${inter_cnt} / ${intra_cnt}))" | bc)
    echo "${graph} ${bench} LLC Misses Reduction by Glign-Inter: ${reduction}"

    echo "========"
done