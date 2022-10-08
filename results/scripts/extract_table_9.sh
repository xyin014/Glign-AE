#!/bin/bash -l

graph=$1
cur_dir=`pwd`
echo $cur_dir

bench=SSSP
for bench in SSSP BFS SSWP SSNP Viterbi Heter; do
    current_bench=$bench
    mode=ligra_s
    file="${cur_dir}/perf/perf_${bench}_${mode}_${graph}.log"
    tmp=`grep -e 'LLC-load-misses:u' $file`
    cnt=`echo ${tmp%'LLC-load-misses:u'*}`
    llcmisses=$(echo "scale=2;(($cnt / 1000000000))" | bc)
    echo "${graph} ${bench} ${mode} LLC Misses: ${llcmisses}"

    mode=ligra_c
    file="${cur_dir}/perf/perf_${bench}_${mode}_${graph}.log"
    tmp=`grep -e 'LLC-load-misses:u' $file`
    cnt=`echo ${tmp%'LLC-load-misses:u'*}`
    llcmisses=$(echo "scale=2;(($cnt / 1000000000))" | bc)
    echo "${graph} ${bench} ${mode} LLC Misses: ${llcmisses}"

    mode=glign
    file="${cur_dir}/perf/perf_${bench}_${mode}_${graph}.log"
    tmp=`grep -e 'LLC-load-misses:u' $file`
    cnt=`echo ${tmp%'LLC-load-misses:u'*}`
    llcmisses=$(echo "scale=2;(($cnt / 1000000000))" | bc)
    echo "${graph} ${bench} ${mode} LLC Misses: ${llcmisses}"

    # mode=glign_intra
    # file="${cur_dir}/perf/perf_${bench}_${mode}_${graph}.log"
    # tmp=`grep -e 'LLC-load-misses:u' $file`
    # cnt=`echo ${tmp%'LLC-load-misses:u'*}`
    # llcmisses=$(echo "scale=2;(($cnt / 1000000000))" | bc)
    # echo "${graph} ${bench} ${mode} LLC Misses: ${llcmisses}"

    # mode=glign_inter
    # file="${cur_dir}/perf/perf_${bench}_${mode}_${graph}.log"
    # tmp=`grep -e 'LLC-load-misses:u' $file`
    # cnt=`echo ${tmp%'LLC-load-misses:u'*}`
    # llcmisses=$(echo "scale=2;(($cnt / 1000000000))" | bc)
    # echo "${graph} ${bench} ${mode} LLC Misses: ${llcmisses}"

    #krill
    if [ "$bench" == "SSSP" ]; then
        bench=M-SSSP
    elif [ "$bench" == "BFS" ]; then
        bench=M-BFS_len
    elif [ "$bench" == "SSWP" ]; then
        bench=M-SSWP
    elif [ "$bench" == "SSNP" ]; then
        bench=M-SSNP
    elif [ "$bench" == "Viterbi" ]; then
        bench=M-Viterbi
    elif [ "$bench" == "Heter" ]; then
        bench=Heter
    fi

    file="${cur_dir}/krill/${graph}/perf_${bench}_krill_${graph}.log"
    tmp=`grep -e 'LLC-load-misses:u' $file`
    cnt=`echo ${tmp%'LLC-load-misses:u'*}`
    llcmisses=$(echo "scale=2;(($cnt / 1000000000))" | bc)
    echo "Krill ${graph} ${bench} LLC Misses: ${llcmisses}"

    #graphm
    # echo "========"
    if [ "$current_bench" == "SSSP" ]; then
        bench=sssp
    elif [ "$current_bench" == "BFS" ]; then
        bench=bfs
    elif [ "$current_bench" == "SSWP" ]; then
        bench=sswp
    elif [ "$current_bench" == "SSNP" ]; then
        bench=ssnp
    elif [ "$current_bench" == "Viterbi" ]; then
        bench=viterbi
    elif [ "$current_bench" == "Heter" ]; then
        bench=heter
    fi

    if [ -e "${cur_dir}/graphm/${graph}/${bench}_perf.txt" ]; then
        `rm "${cur_dir}/graphm/${graph}/${bench}_perf.txt"`
    fi
    cat ${cur_dir}/graphm/${graph}/perf_concurrent_${bench}_graphm*.log >> "${cur_dir}/graphm/${graph}/${bench}_perf.txt"
    file="${cur_dir}/graphm/${graph}/${bench}_perf.txt"
    graphm_cnt=`grep -e 'LLC-load-misses:u' $file | awk '{SUM += $(1)} END {print SUM }'`
    # echo $graphm_cnt
    llcmisses=$(echo "scale=2;(($graphm_cnt / 1000000000))" | bc)
    echo "GraphM ${graph} ${bench} LLC Misses: ${llcmisses}"

    echo "================"
done