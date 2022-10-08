#!/bin/bash -l
date

# if on HPCC
# module load extra/1
# module load GCCcore/7.4.0

name=$1
graph_path=$2
query_file=$3
output_path=../results/${name}

# name=LJ
# graph_path=~/bigdata/data/soc-LiveJournal1.weighted.adj
# query_file=~/bigdata/query_inputs/lj_hops_evenly_shuffled.txt
# output_path=~/bigdata/asplos23-ae/LJ

# graph_path=~/bigdata/data/twitter-2010.weighted.adj
# query_file=~/bigdata/query_inputs/twitter_hops_evenly_512.txt
# output_path=~/bigdata/asplos23-ae

total=512
# SSSP
bench=SSSP
# running ligra-s
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 1 -qf $query_file $graph_path > ${output_path}/${bench}_ligra_s_${name}.txt
# running glign
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_${name}.txt
# running ligra-c
./${bench}_Batch -option ligra-c -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_ligra_c_${name}.txt
# running glign-intra
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -qf $query_file $graph_path > ${output_path}/${bench}_glign_intra_${name}.txt
# running glign-inter
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_inter_${name}.txt
# running glign-batch
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -qf $query_file $graph_path > ${output_path}/${bench}_glign_batch_${name}.txt

# BFS
bench=BFS
# running ligra-s
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 1 -qf $query_file $graph_path > ${output_path}/${bench}_ligra_s_${name}.txt
# running glign
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_${name}.txt
# running ligra-c
./${bench}_Batch -option ligra-c -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_ligra_c_${name}.txt
# running glign-intra
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -qf $query_file $graph_path > ${output_path}/${bench}_glign_intra_${name}.txt
# running glign-inter
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_inter_${name}.txt
# running glign-batch
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -qf $query_file $graph_path > ${output_path}/${bench}_glign_batch_${name}.txt

# SSWP
bench=SSWP
# running ligra-s
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 1 -qf $query_file $graph_path > ${output_path}/${bench}_ligra_s_${name}.txt
# running glign
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_${name}.txt
# running ligra-c
./${bench}_Batch -option ligra-c -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_ligra_c_${name}.txt
# running glign-intra
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -qf $query_file $graph_path > ${output_path}/${bench}_glign_intra_${name}.txt
# running glign-inter
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_inter_${name}.txt
# running glign-batch
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -qf $query_file $graph_path > ${output_path}/${bench}_glign_batch_${name}.txt

# SSNP
bench=SSNP
# running ligra-s
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 1 -qf $query_file $graph_path > ${output_path}/${bench}_ligra_s_${name}.txt
# running glign
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_${name}.txt
# running ligra-c
./${bench}_Batch -option ligra-c -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_ligra_c_${name}.txt
# running glign-intra
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -qf $query_file $graph_path > ${output_path}/${bench}_glign_intra_${name}.txt
# running glign-inter
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_inter_${name}.txt
# running glign-batch
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -qf $query_file $graph_path > ${output_path}/${bench}_glign_batch_${name}.txt

# Viterbi
bench=Viterbi
# running ligra-s
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 1 -qf $query_file $graph_path > ${output_path}/${bench}_ligra_s_${name}.txt
# running glign
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_${name}.txt
# running ligra-c
./${bench}_Batch -option ligra-c -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_ligra_c_${name}.txt
# running glign-intra
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -qf $query_file $graph_path > ${output_path}/${bench}_glign_intra_${name}.txt
# running glign-inter
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 2 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_inter_${name}.txt
# running glign-batch
./${bench}_Batch -option glign -batch 64 -max_combination $total -mode 3 -qf $query_file $graph_path > ${output_path}/${bench}_glign_batch_${name}.txt

# Heter
bench=Heter
# running ligra-s
./${bench}_Batch -option glign-heter -batch 64 -max_combination $total -mode 1 -qf $query_file $graph_path > ${output_path}/${bench}_ligra_s_${name}.txt
# running glign
./${bench}_Batch -option glign-heter -batch 64 -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_${name}.txt
# running ligra-c
./${bench}_Batch -option ligra-c -batch 64 -max_combination $total -qf $query_file $graph_path > ${output_path}/${bench}_ligra_c_${name}.txt
# running glign-intra
./${bench}_Batch -option glign-heter -batch 64 -max_combination $total -mode 2 -qf $query_file $graph_path > ${output_path}/${bench}_glign_intra_${name}.txt
# running glign-inter
./${bench}_Batch -option glign-heter -batch 64 -max_combination $total -mode 2 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_inter_${name}.txt
# running glign-batch
./${bench}_Batch -option glign-heter -batch 64 -max_combination $total -mode 3 -qf $query_file $graph_path > ${output_path}/${bench}_glign_batch_${name}.txt


# iBFS
bench=BFS
./${bench}_Batch -option ibfs -batch 64 -max_combination $total -qf $query_file $graph_path > ${output_path}/${bench}_ibfs_${name}.txt


# Varying batch sizes
# SSSP
total=256
bench=SSSP
for bsize in 2 4 8 16 32 64 128; do
    echo $bsize
    # running glign-batch
    ./${bench}_Batch -option glign -batch $bsize -max_combination $total -mode 3 -qf $query_file $graph_path > ${output_path}/${bench}_glign_batch_${name}_${bsize}.txt
    # running gling-inter
    ./${bench}_Batch -option glign -batch $bsize -max_combination $total -mode 2 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_inter_${name}_${bsize}.txt
    # running ligra-c
    ./${bench}_Batch -option ligra-c -batch $bsize -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_ligra_c_${name}_${bsize}.txt
done

bench=BFS
for bsize in 2 4 8 16 32 64 128; do
    echo $bsize
    # running glign-batch
    ./${bench}_Batch -option glign -batch $bsize -max_combination $total -mode 3 -qf $query_file $graph_path > ${output_path}/${bench}_glign_batch_${name}_${bsize}.txt
    # running gling-inter
    ./${bench}_Batch -option glign -batch $bsize -max_combination $total -mode 2 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_inter_${name}_${bsize}.txt
    # running ligra-c
    ./${bench}_Batch -option ligra-c -batch $bsize -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_ligra_c_${name}_${bsize}.txt
done 

bench=SSWP
for bsize in 2 4 8 16 32 64 128; do
    echo $bsize
    # running glign-batch
    ./${bench}_Batch -option glign -batch $bsize -max_combination $total -mode 3 -qf $query_file $graph_path > ${output_path}/${bench}_glign_batch_${name}_${bsize}.txt
    # running gling-inter
    ./${bench}_Batch -option glign -batch $bsize -max_combination $total -mode 2 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_inter_${name}_${bsize}.txt
    # running ligra-c
    ./${bench}_Batch -option ligra-c -batch $bsize -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_ligra_c_${name}_${bsize}.txt
done

bench=Viterbi
for bsize in 2 4 8 16 32 64 128; do
    echo $bsize
    # running glign-batch
    ./${bench}_Batch -option glign -batch $bsize -max_combination $total -mode 3 -qf $query_file $graph_path > ${output_path}/${bench}_glign_batch_${name}_${bsize}.txt
    # running gling-inter
    ./${bench}_Batch -option glign -batch $bsize -max_combination $total -mode 2 -delay -qf $query_file $graph_path > ${output_path}/${bench}_glign_inter_${name}_${bsize}.txt
    # running ligra-c
    ./${bench}_Batch -option ligra-c -batch $bsize -max_combination $total -mode 3 -delay -qf $query_file $graph_path > ${output_path}/${bench}_ligra_c_${name}_${bsize}.txt
done

date
