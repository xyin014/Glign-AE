#!/bin/bash -l

graph=$1
cur_dir=`pwd`

bench=SSSP
for bench in SSSP BFS SSWP SSNP Viterbi Heter; do
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

    mode=glign
    file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
    glign=`grep -e 'Glign evaluation time: ' $file | awk '{print $NF}'`
    speedup=$(echo "scale=2;(($seq / $glign))" | bc)
    echo "${graph} ${bench} ${mode} speedup: ${speedup}"

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
    file="${cur_dir}/krill/${graph}/${bench}_${graph}.txt"
    krill_time=`grep -e 'Running time: ' $file | awk '{print $NF}' | grep -Eo '[+-]?[0-9]+([.][0-9]+)?'`
    speedup=$(echo "scale=2;(($seq / ${krill_time}))" | bc)
    echo "Krill ${graph} ${bench} speedup: ${speedup}"

    #graphm
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

    if [ -e "${cur_dir}/graphm/${graph}/${bench}_res.txt" ]; then
        `rm "${cur_dir}/graphm/${graph}/${bench}_res.txt"`
    fi
    cat ${cur_dir}/graphm/${graph}/concurrent_${bench}_graphm*.txt >> "${cur_dir}/graphm/${graph}/${bench}_res.txt"
    file="${cur_dir}/graphm/${graph}/${bench}_res.txt"
    graphm_time=`grep -e 'iterations of concurrent jobs took ' $file | awk '{SUM += $(NF-1)} END {print SUM }'`
    # {print $(NF-1)}'`
    # echo $graphm_time
    speedup=$(echo "scale=2;(($seq / ${graphm_time}))" | bc)
    echo "GraphM ${graph} ${bench} speedup: ${speedup}"

    echo "========"
done