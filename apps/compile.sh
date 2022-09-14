#!/bin/bash -l

# git checkout ae-checkin
make clean
CILK=1 EDGELONG=1 LONG=1 make SSSP_Batch
CILK=1 EDGELONG=1 LONG=1 make BFS_Batch
CILK=1 EDGELONG=1 LONG=1 make SSWP_Batch
CILK=1 EDGELONG=1 LONG=1 make Viterbi_Batch
CILK=1 EDGELONG=1 LONG=1 make SSNP_Batch

# git checkout ae-heter
# CILK=1 EDGELONG=1 LONG=1 make Heter_Batch
# git checkout ae-checkin
