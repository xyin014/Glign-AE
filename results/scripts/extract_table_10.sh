#!/bin/bash -l

graph=$1
cur_dir=`pwd`

bench=SSSP
for bench in SSSP BFS SSWP SSNP Viterbi Heter; do
    
    mode=ligra_c
    file="${cur_dir}/perf/perf_${bench}_${mode}_${graph}.log"
    tmp=`grep -e 'LLC-load-misses:u' $file`
    ligra_c_cnt=`echo ${tmp%'LLC-load-misses:u'*}`
    # llcmisses=$(echo "scale=2;(($ligra_c_cnt / 1000000000))" | bc)
    # echo "${graph} ${bench} ${mode} LLC Misses: ${llcmisses}"

    mode=glign_intra
    file="${cur_dir}/perf/perf_${bench}_${mode}_${graph}.log"
    tmp=`grep -e 'LLC-load-misses:u' $file`
    intra_cnt=`echo ${tmp%'LLC-load-misses:u'*}`
    reduction=$(echo "scale=2;((${intra_cnt} / ${ligra_c_cnt}))" | bc)
    echo "${graph} ${bench} LLC Misses Reduction by Glign-Intra: ${reduction}"

    echo "========"
done