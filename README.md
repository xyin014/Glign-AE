Glign: Taming Misaligned Graph Traversals in Concurrent Graph Processing
======================

This repository provides experiments and scripts for reproducing most figures and tables in our paper *Glign: Taming Misaligned Graph Traversals in Concurrent Graph Processing*. Our paper introduces *Glign*, a runtime system for evaluating concurrent graph queries. Currently, Glign is built on top of [Ligra](https://github.com/jshun/ligra). 

In the instructions below, we provide steps for compiling and running experiments.

Compilation
--------

Compilation is done from within the apps/ directory. The compiled code
will work on both uncompressed and compressed graphs and hypergraphs.

Compilers

* g++ &gt;= 5.3.0 and &lt;= 7.4.0 with support for Cilk Plus

After the appropriate environment variables are set, to compile,
first checkout `ae-checkin` branch and run

```
$ bash compile.sh
```
This will compile all six benchmakrs in the paper (SSSP_Batch, BFS_Batch, SSWP_Batch, SSNP_Batch, Viterbi_Batch, Heter_Batch).

Running concurrent queries in Glign.
-------
For example, to run concurrent SSSP queries:

```
$ ./SSSP_Batch -option glign -batch 64 -max_combination 512 -mode 3 -delay -qf [queries source file] [graph]
``` 
The above will evaluate 512 queries in 8 batches, each batch contains 64 concurrent queries.
Query files are provided (in `query_input` folder), the queries for each graph are generated based on Section 4.1 of the paper.
`-batch` denotes the batch size, `-max_combination` denotes the total number of queries to evaluate. [graph] should be weighted graph in adjacency format. User can use Ligra's tool to convert SNAP graph format to adjacency format and add random weights.

There are some variants of Glign in the paper which are summarized below:
| **Variant** | **parameters**               |
|-------------|------------------------------|
| Ligra-Seq   | -option glign -mode 1        |
| Ligra-C     | -option ligra-c              |
| Glign-Intra | -option glign -mode 2        |
| Glign-Inter | -option glign -mode 2 -delay |
| Glign-Batch | -option glign -mode 3        |
| Glign       | -option glign -mode 3 -delay |

For heterogeneous queries (a batch of mixed SSSP, SSWP, BFS, and SSNP):
| **Variant** | **parameters**                     |
|-------------|------------------------------------|
| Ligra-Seq   | -option glign-heter -mode 1        |
| Ligra-C     | -option ligra-c                    |
| Glign-Intra | -option glign-heter -mode 2        |
| Glign-Inter | -option glign-heter -mode 2 -delay |
| Glign-Batch | -option glign-heter -mode 3        |
| Glign       | -option glign-heter -mode 3 -delay |

In the future release, we will simpilify the config parameters for different settings.

Reproducing figures and tables in the paper
-------
After converting SNAP graphs to adjacency graphs for Ligra and Krill, run the following commands (under the `apps` folder):

The overall performance tests in Glign:
```
bash run_experiments.sh [name] [graph_path] [query_file]
```
[name] should be the abbreviation of tested graphs in the paper: LJ, WP, UK2, TW, FR, RDCA, RDUS.
[graph_path] is the location of weighted adjacency graph.
[query_file] is the query input (in `../query_input/`).
The results are stored in `../results/[name]`.

To run profiling experiments in Glign:
```
bash run_profilings.sh [name] [graph_path] [query_file]
```

To run performance tests in Krill (under `Krill-AE/apps` directory):
```
bash compile_benchmarks.sh
bash run_experiments.sh [name] [graph_path] [query_file] [output_path]
```
Note that [output_path] is the root directory of Glign (`.../Glign-AE`).

To run profiling experiments in **Krill** (under `Krill-AE/apps` directory):
```
bash run_profilings.sh [name] [graph_path] [query_file] [output_path]
```

To run performance tests in **GraphM** (under `GraphM-AE` directory):
```
bash run_experiments.sh [name] [graph_path] [query_file] [output_path] [graph_size]
```
Note that [graph_size] is a required input parameter by GraphM. Use 828 for LJ, 5247 for WP, 17621 for TW, 21673 for FR, 3142 for UK2.
[output_path] is the root directory of Glign (`.../Glign-AE`).

To run profiling experiments in GraphM (under `GraphM-AE` directory):
```
bash run_profilings.sh [name] [graph_path] [query_file] [output_path] [graph_size]
```
*Generating data for figures and tables*
Under the `.../Glign-AE/results` directory and run the following command to:

Generate data for Figure 11:
```
bash ./scripts/extract_figure_11.sh [graph_abbr]
```
[graph_abbr] is the name of graph, use LJ for LiveJournal.

Generate data for Figure 12:
```
bash ./scripts/extract_figure_12.sh [graph_abbr]
```

Generate data for Figure 13:
```
bash ./scripts/extract_figure_13.sh [graph_abbr]
```

Generate data for Figure 14:
```
bash ./scripts/extract_figure_14.sh [graph_abbr]
```

Generate data for Figure 15:
```
bash ./scripts/extract_figure_15.sh [graph_abbr]
```

Generate data for Figure 16:
```
bash ./scripts/extract_figure_16.sh [graph_abbr]
```

Generate data for Table 8:
```
bash ./scripts/extract_table_8.sh [graph_abbr]
```

Generate data for Table 9:
```
bash ./scripts/extract_table_9.sh [graph_abbr]
```
The cache misses are collected by `perf` tool, please make sure the CPU supports related hardware counters.

Generate data for Table 10:
```
bash ./scripts/extract_table_10.sh [graph_abbr]
```

Generate data for Table 12:
```
bash ./scripts/extract_table_12.sh [graph_abbr]
```

Generate data for Table 13:
```
../apps/SSSP_Batch -option ground-truth -batch 2 -max_combination 512 -qf ../query_input/LJ_queries.txt ../ae-data/soc-LiveJournal1.weighted.adj
```
The results table will be printed in the terminal after the program finishes.

Data for Table 14:
Please check any of Glign's output file inside the results folder, there is a line `Profiling cost: ` in each file.

Data for Table 15:
This is similar to Figure 11.
