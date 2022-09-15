#!/bin/bash -l

# dir=$1
graph=$1
mode=ligra_s
# for graph in LJ WP UK2 TW FR; do
bench=SSSP
cur_dir=`pwd`
file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
# echo "file: ${file}"
sssp_lj_seq=`grep -e 'Sequential evaluation time: ' $file | awk '{print $NF}'`
echo "${graph} ${bench} ${mode} seq time: ${sssp_lj_seq}"

bench=BFS
cur_dir=`pwd`
file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
# echo "file: ${file}"
bfs_lj_seq=`grep -e 'Sequential evaluation time: ' $file | awk '{print $NF}'`
echo "${graph} ${bench} ${mode} seq time: ${bfs_lj_seq}"

bench=SSWP
cur_dir=`pwd`
file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
# echo "file: ${file}"
sswp_lj_seq=`grep -e 'Sequential evaluation time: ' $file | awk '{print $NF}'`
echo "${graph} ${bench} ${mode} seq time: ${sswp_lj_seq}"

bench=SSNP
cur_dir=`pwd`
file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
# echo "file: ${file}"
ssnp_lj_seq=`grep -e 'Sequential evaluation time: ' $file | awk '{print $NF}'`
echo "${graph} ${bench} ${mode} seq time: ${ssnp_lj_seq}"

bench=Viterbi
cur_dir=`pwd`
file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
# echo "file: ${file}"
viterbi_lj_seq=`grep -e 'Sequential evaluation time: ' $file | awk '{print $NF}'`
echo "${graph} ${bench} ${mode} seq time: ${viterbi_lj_seq}"

bench=Heter
cur_dir=`pwd`
file="${cur_dir}/${graph}/${bench}_${mode}_${graph}.txt"
# echo "file: ${file}"
heter_lj_seq=`grep -e 'Sequential evaluation time: ' $file | awk '{print $NF}'`
echo "${graph} ${bench} ${mode} seq time: ${heter_lj_seq}"
# done