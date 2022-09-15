#!/bin/bash -l
date

# if on HPCC
# module load extra/1
# module load GCCcore/7.4.0

name=$1
graph_path=$2
query_file=$3
output_path=../results/perf

echo $name
# name=LJ
# graph_path=~/bigdata/data/soc-LiveJournal1.weighted.adj
# query_file=~/bigdata/query_inputs/lj_hops_evenly_shuffled.txt
# output_path=~/bigdata/asplos23-ae/perf

total=512

# SSSP
bench=SSSP
# running ligra-s
perf stat -d -d -o $output_path/"perf_${bench}_ligra_s_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 1 -qf $query_file $graph_path
# running glign
perf stat -d -d -o $output_path/"perf_${bench}_glign_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path
# running ligra-c
perf stat -d -d -o $output_path/"perf_${bench}_ligra_c_${name}.log" ./${bench}_Batch -option ligra-c -batch 64 -max_combination $total -qf $query_file $graph_path
# running glign-intra
perf stat -d -d -o $output_path/"perf_${bench}_glign_intra_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -qf $query_file $graph_path
# running glign-inter
perf stat -d -d -o $output_path/"perf_${bench}_glign_inter_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -delay -qf $query_file $graph_path

bench=BFS
# running ligra-s
perf stat -d -d -o $output_path/"perf_${bench}_ligra_s_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 1 -qf $query_file $graph_path
# running glign
perf stat -d -d -o $output_path/"perf_${bench}_glign_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path
# running ligra-c
perf stat -d -d -o $output_path/"perf_${bench}_ligra_c_${name}.log" ./${bench}_Batch -option ligra-c -batch 64 -max_combination $total -qf $query_file $graph_path
# running glign-intra
perf stat -d -d -o $output_path/"perf_${bench}_glign_intra_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -qf $query_file $graph_path
# running glign-inter
perf stat -d -d -o $output_path/"perf_${bench}_glign_inter_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -delay -qf $query_file $graph_path

bench=SSWP
# running ligra-s
perf stat -d -d -o $output_path/"perf_${bench}_ligra_s_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 1 -qf $query_file $graph_path
# running glign
perf stat -d -d -o $output_path/"perf_${bench}_glign_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path
# running ligra-c
perf stat -d -d -o $output_path/"perf_${bench}_ligra_c_${name}.log" ./${bench}_Batch -option ligra-c -batch 64 -max_combination $total -qf $query_file $graph_path
# running glign-intra
perf stat -d -d -o $output_path/"perf_${bench}_glign_intra_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -qf $query_file $graph_path
# running glign-inter
perf stat -d -d -o $output_path/"perf_${bench}_glign_inter_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -delay -qf $query_file $graph_path

bench=SSNP
# running ligra-s
perf stat -d -d -o $output_path/"perf_${bench}_ligra_s_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 1 -qf $query_file $graph_path
# running glign
perf stat -d -d -o $output_path/"perf_${bench}_glign_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path
# running ligra-c
perf stat -d -d -o $output_path/"perf_${bench}_ligra_c_${name}.log" ./${bench}_Batch -option ligra-c -batch 64 -max_combination $total -qf $query_file $graph_path
# running glign-intra
perf stat -d -d -o $output_path/"perf_${bench}_glign_intra_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -qf $query_file $graph_path
# running glign-inter
perf stat -d -d -o $output_path/"perf_${bench}_glign_inter_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -delay -qf $query_file $graph_path

bench=Viterbi
# running ligra-s
perf stat -d -d -o $output_path/"perf_${bench}_ligra_s_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 1 -qf $query_file $graph_path 
# running glign
perf stat -d -d -o $output_path/"perf_${bench}_glign_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path 
# running ligra-c
perf stat -d -d -o $output_path/"perf_${bench}_ligra_c_${name}.log" ./${bench}_Batch -option ligra-c -batch 64 -max_combination $total -qf $query_file $graph_path 
# running glign-intra
perf stat -d -d -o $output_path/"perf_${bench}_glign_intra_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -qf $query_file $graph_path 
# running glign-inter
perf stat -d -d -o $output_path/"perf_${bench}_glign_inter_${name}.log" ./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -delay -qf $query_file $graph_path 


# Heter
bench=Heter
# running ligra-s
perf stat -d -d -o $output_path/"perf_${bench}_ligra_s_${name}.log" ./${bench}_Batch -option glign-heter -batch 64 -max_combination $total -mode 1 -qf $query_file $graph_path 
# running glign
perf stat -d -d -o $output_path/"perf_${bench}_glign_${name}.log" ./${bench}_Batch -option glign-heter -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path
# running ligra-c
perf stat -d -d -o $output_path/"perf_${bench}_ligra_c_${name}.log" ./${bench}_Batch -option ligra-c -batch 64 -max_combination $total -qf $query_file $graph_path 
# running glign-intra
perf stat -d -d -o $output_path/"perf_${bench}_glign_intra_${name}.log" ./${bench}_Batch -option glign-heter -batch 64 -max_combination $total -mode 2 -qf $query_file $graph_path 
# running glign-inter
perf stat -d -d -o $output_path/"perf_${bench}_glign_intra_${name}.log" ./${bench}_Batch -option glign-heter -batch 64 -max_combination $total -mode 2 -delay -qf $query_file $graph_path