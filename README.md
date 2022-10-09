# Glign: Taming Misaligned Graph Traversals in Concurrent Graph Processing
======================

This repository provides experiments and scripts for reproducing most figures and tables in our paper *Glign: Taming Misaligned Graph Traversals in Concurrent Graph Processing*. Our paper introduces *Glign*, a runtime system for evaluating concurrent graph queries. Currently, Glign is built on top of [Ligra](https://github.com/jshun/ligra). 

In the instructions below, we provide steps for compiling and running experiments.

## Dependency requirements
--------
G++ compilers with support for Cilk Plus is required (version &gt;= 5.3.0 and &lt;= 7.4.0).
[Perf tool](https://man7.org/linux/man-pages/man1/perf.1.html) is used for profiling.

## Compilation
--------
Compilers

* g++ &gt;= 5.3.0 and &lt;= 7.4.0 with support for Cilk Plus

To compile, first checkout `ae-checkin` branch, cd into `Glign-AE/apps` folder and run

```
$ bash compile.sh
```
This will compile all six benchmakrs in the paper (SSSP_Batch, BFS_Batch, SSWP_Batch, SSNP_Batch, Viterbi_Batch, Heter_Batch).

## Glign usages.
-------
For example, to run concurrent SSSP queries (under `Glign-AE/apps` directory):

```
$ ./SSSP_Batch -option glign -batch 64 -max_combination 512 -mode 3 -delay -qf [queries source file] [graph]
``` 
The above will evaluate 512 queries in 8 batches, each batch contains 64 concurrent queries.
Query files are provided (in `Glign-AE/query_input` folder), the queries for each graph are generated based on Section 4.1 of the paper.
`-batch` denotes the batch size, `-max_combination` denotes the total number of queries to evaluate. [graph] should be a weighted graph in adjacency format. User can use Ligra's tool to convert SNAP graph format to adjacency format and add random weights.
A smaller graph (LiveJournal) has been provided inside the `ae-data` folder. The GridGrid graph format (used by GraphM) for LJ is also provided in `ae-data/graphm`.

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

## Running experiments
-------
The adj graphs (for Ligra and Krill) and grid format graphs (for GraphM) can be downloaded from the [shared Google Drive folder](https://drive.google.com/drive/folders/1VNzred7_cdvoyvwyQdBnQkLJZqaDFYuV?usp=sharing).

**It is suggested that first finish data collection for one graph. If time allows, larger graphs can be further tested.**

| **Adj Graph** | **File name**               | **Size** |
|-------------|------------------------------|-----------|
| LJ   | soc-LiveJournal1.weighted.adj       | 0.69GB |
| WP     | wikipedia_link_en_weighted.adj    | 4.3GB  |
| UK2 | uk-2002-sym-weighted.adj        | 5.5 GB |
| TW | twitter-2010.weighted.adj | 16.1GB |
| FR | friendster.sym.weighted.adj  | 40.3GB |
| RD-CA       | roadNet-CA.weighted.adj | 0.067GB |
| RD-US       | road-road-usa.weighted.adj | 0.82GB |

<!-- After converting SNAP graphs to adjacency graphs for Ligra and Krill,  -->
Note: it is possible to download the SNAP graphs (from [SNAP Project](https://snap.stanford.edu/data/index.html), and convert SNAP graphs to Adj format (following notes in https://github.com/jshun/ligra#graph-converters), and convert SNAP graphs to Grid format (following notes in https://github.com/Program-Anonymity/GraphM)), 
however, the Goolge Drive shared folder contains all necessary graph inputs. 


After downloading the graph files, please copy and paste all files into `.../Glign/ae-data` folder and
run the following commands (under the `.../Glign/apps` folder):

**1. The overall performance tests in Glign (in the `.../Glign/apps` folder):**
```
bash run_experiments.sh [name] [graph_path] [query_file]
```
[name] should be the abbreviation of tested graphs in the paper: LJ, WP, UK2, TW, FR, RDCA, RDUS.
[graph_path] is the location of weighted adjacency graph.
[query_file] is the query input (`../query_input/[name]`).
The results are stored in `../results/[name]`.

Example commands for running different graphs are listed below:
<!-- bash run_experiments.sh LJ ../ae-data/soc-LiveJournal1.weighted.adj ../query_input/LJ_queries.txt -->
| **Graph** | **Command** | **Estimated time** |
|-------------|------------------------------|-------------|
| LJ    | bash run_experiments.sh LJ ../ae-data/soc-LiveJournal1.weighted.adj ../query_input/LJ_queries.txt | less than 2 hours|
| WP    | bash run_experiments.sh WP ../ae-data/wikipedia_link_en_weighted.adj ../query_input/WP_queries.txt | ~4 hours|
| UK2   | bash run_experiments.sh UK2 ../ae-data/uk-2002-sym-weighted.adj ../query_input/UK2_queries.txt | ~5 hours|
| TW    | bash run_experiments.sh TW ../ae-data/twitter-2010.weighted.adj ../query_input/TW_queries.txt | ~12 hours|
| FR    | bash run_experiments.sh FR ../ae-data/friendster.sym.weighted.adj ../query_input/FR_queries.txt | > 24 hours|
<!-- | RD-CA | bash run_experiments.sh RDCA ../ae-data/roadNet-CA.weighted.adj ../query_input/road-ca-512-random.txt | ~2 hours| -->
<!-- | RD-US | bash run_experiments.sh RDUS ../ae-data/road-road-usa.weighted.adj ../query_input/road-usa-512-random.txt | ~8 hours| -->

**2. To run road networks (in the `.../Glign/apps` folder) (only for Table 15):**
| **Graph** | **Command** | **Estimated time** |
|-------------|------------------------------|-------------|
| RD-CA | bash run_roadnets.sh RDCA ../ae-data/roadNet-CA.weighted.adj ../query_input/road-ca-512-random.txt | ~2 hours|
| RD-US | bash run_roadnets.sh RDUS ../ae-data/road-road-usa.weighted.adj ../query_input/road-usa-512-random.txt | ~8 hours|


**3. To run profiling experiments in Glign:**
```
bash run_profilings.sh [name] [graph_path] [query_file]
```
For example, to run experiments for LJ, use `bash run_profilings.sh LJ ../ae-data/soc-LiveJournal1.weighted.adj ../query_input/LJ_queries.txt`

Example commands for running different graphs are listed below:
| **Graph** | **Command** | **Estimated time** |
|-------------|------------------------------|-------------|
| LJ    | bash run_profilings.sh LJ ../ae-data/soc-LiveJournal1.weighted.adj ../query_input/LJ_queries.txt | less than 2 hours|
| WP    | bash run_profilings.sh WP ../ae-data/wikipedia_link_en_weighted.adj ../query_input/WP_queries.txt |~4 hours|
| UK2   | bash run_profilings.sh UK2 ../ae-data/uk-2002-sym-weighted.adj ../query_input/UK2_queries.txt | ~5 hours|
| TW    | bash run_profilings.sh TW ../ae-data/twitter-2010.weighted.adj ../query_input/TW_queries.txt | ~12 hours|
| FR    | bash run_profilings.sh FR ../ae-data/friendster.sym.weighted.adj ../query_input/FR_queries.txt | >24 hours|
<!-- | RD-CA | bash run_profilings.sh RDCA ../ae-data/roadNet-CA.weighted.adj ../query_input/road-ca-512-random.txt | ~2 hours| -->
<!-- | RD-US | bash run_profilings.sh RDUS ../ae-data/road-road-usa.weighted.adj ../query_input/road-usa-512-random.txt | ~8 hours| -->

**4. To run performance tests in Krill (under `Krill-AE/apps` directory):**
```
bash compile_benchmarks.sh
bash run_experiments.sh [name] [graph_path] [query_file] [output_path]
```
Note that [output_path] is the root directory of Glign (`.../Glign-AE`).

Example commands for running different graphs are listed below:
| **Graph** | **Command** | **Estimated time** |
|-------------|------------------------------|-------------|
| LJ    | bash run_experiments.sh LJ ../../Glign-AE/ae-data/soc-LiveJournal1.weighted.adj ../../Glign-AE/query_input/LJ_queries.txt ../../Glign-AE | less than 2 hours|
| WP    | bash run_experiments.sh WP ../../Glign-AE/ae-data/wikipedia_link_en_weighted.adj ../../Glign-AE/query_input/WP_queries.txt ../../Glign-AE | ~4 hours|
| UK2   | bash run_experiments.sh UK2 ../../Glign-AE/ae-data/uk-2002-sym-weighted.adj ../../Glign-AE/query_input/UK2_queries.txt ../../Glign-AE | ~6 hours|
| TW    | bash run_experiments.sh TW ../../Glign-AE/ae-data/twitter-2010.weighted.adj ../../Glign-AE/query_input/TW_queries.txt ../../Glign-AE | ~12 hours|
| FR    | bash run_experiments.sh FR ../../Glign-AE/ae-data/friendster.sym.weighted.adj ../../Glign-AE/query_input/FR_queries.txt ../../Glign-AE | >24 hours|

**5. To run profiling experiments in **Krill** (under `Krill-AE/apps` directory):**
```
bash run_profilings.sh [name] [graph_path] [query_file] [output_path]
```

Example commands for running different graphs are listed below:
| **Graph** | **Command** | **Estimated time** |
|-------------|------------------------------|-------------|
| LJ    | bash run_profilings.sh LJ ../../Glign-AE/ae-data/soc-LiveJournal1.weighted.adj ../../Glign-AE/query_input/LJ_queries.txt ../../Glign-AE | less than 2 hours|
| WP    | bash run_profilings.sh WP ../../Glign-AE/ae-data/wikipedia_link_en_weighted.adj ../../Glign-AE/query_input/WP_queries.txt ../../Glign-AE | ~4 hours|
| UK2   | bash run_profilings.sh UK2 ../../Glign-AE/ae-data/uk-2002-sym-weighted.adj ../../Glign-AE/query_input/query_input/UK2_queries.txt ../../Glign-AE | ~6 hours|
| TW    | bash run_profilings.sh TW ../../Glign-AE/ae-data/twitter-2010.weighted.adj ../../Glign-AE/query_input/TW_queries.txt ../../Glign-AE | ~12 hours|
| FR    | bash run_profilings.sh FR ../../Glign-AE/ae-data/friendster.sym.weighted.adj ../../Glign-AE/query_input/FR_queries.txt ../../Glign-AE | >24 hours|


**6. To run performance tests in **GraphM** (under `GraphM-AE` directory):**
```
bash compile_benchmarks.sh
bash run_experiments.sh [name] [graph_path] [query_file] [output_path] [graph_size]
```
Note that [graph_size] is a required input parameter by GraphM. Use 828 for LJ, 5247 for WP, 17621 for TW, 21673 for FR, 3142 for UK2.
[output_path] is the root directory of Glign (`../Glign-AE`).

Example commands for running different graphs are listed below:
| **Graph** | **Command** | **Estimated time** |
|-------------|------------------------------|-------------|
| LJ    | bash run_experiments.sh LJ ../Glign-AE/ae-data/graphm/LJ ../Glign-AE/query_input/LJ_queries.txt ../Glign-AE 828 | less than 2 hours|
| WP    | bash run_experiments.sh WP ../Glign-AE/ae-data/graphm/WP ../Glign-AE/query_input/WP_queries.txt ../Glign-AE 5247 | ~4 hours|
| UK2   | bash run_experiments.sh UK2 ../Glign-AE/ae-data/graphm/UK2 ../Glign-AE/query_input/UK2_queries.txt ../Glign-AE 3142 | ~6 hours|
| TW    | bash run_experiments.sh TW ../Glign-AE/ae-data/graphm/TW ../Glign-AE/query_input/TW_queries.txt ../Glign-AE 17621 | >12 hours|
| FR    | bash run_experiments.sh FR ../Glign-AE/ae-data/graphm/FR ../Glign-AE/query_input/FR_queries.txt ../Glign-AE 21673 | >24 hours|

**7. To run profiling experiments in GraphM (under `GraphM-AE` directory):**
```
bash run_profilings.sh [name] [graph_path] [query_file] [output_path] [graph_size]
```

Example commands for running different graphs are listed below:
| **Graph** | **Command** | **Estimated time** |
|-------------|------------------------------|-------------|
| LJ    | bash run_profilings.sh LJ ../Glign-AE/ae-data/graphm/LJ ../Glign-AE/query_input/LJ_queries.txt ../Glign-AE 828 | less than 2 hours|
| WP    | bash run_profilings.sh WP ../Glign-AE/ae-data/graphm/WP ../Glign-AE/query_input/WP_queries.txt ../Glign-AE 5247 | ~4 hours|
| UK2   | bash run_profilings.sh UK2 ../Glign-AE/ae-data/graphm/UK2 ../Glign-AE/query_input/UK2_queries.txt ../Glign-AE 3142 | ~6 hours|
| TW    | bash run_profilings.sh TW ../Glign-AE/ae-data/graphm/TW ../Glign-AE/query_input/TW_queries.txt ../Glign-AE 17621 | >12 hours|
| FR    | bash run_profilings.sh FR ../Glign-AE/ae-data/graphm/FR ../Glign-AE/query_input/FR_queries.txt ../Glign-AE 21673 | >24 hours|


## Generating data for figures and tables

Under the `.../Glign-AE/results` directory and run the following command to:

**1. Generate data for Figure 11:**
```
bash ./scripts/extract_figure_11.sh [graph_abbr]
```
[graph_abbr] is the name of graph, use LJ for LiveJournal.
| **graph_abbr** | **Command** |
|-------------|------------------------------|
| LJ    | bash ./scripts/extract_figure_11.sh LJ |
| WP    | bash ./scripts/extract_figure_11.sh WP |
| UK2   | bash ./scripts/extract_figure_11.sh UK2 |
| TW    | bash ./scripts/extract_figure_11.sh TW |
| FR    | bash ./scripts/extract_figure_11.sh FR |
<!-- | RD-CA | bash ./scripts/extract_figure_11.sh RDCA |
| RD-US | bash ./scripts/extract_figure_11.sh RDUS | -->

**2. Generate data for Figure 12:**
```
bash ./scripts/extract_figure_12.sh [graph_abbr]
```

**3. Generate data for Figure 13:**
```
bash ./scripts/extract_figure_13.sh [graph_abbr]
```

**4. Generate data for Figure 14:**
```
bash ./scripts/extract_figure_14.sh [graph_abbr]
```

**5. Generate data for Figure 15:**
```
bash ./scripts/extract_figure_15.sh [graph_abbr]
```

**6. Generate data for Figure 16:**
```
bash ./scripts/extract_figure_16.sh [graph_abbr]
```

**7. Generate data for Table 8:**
```
bash ./scripts/extract_table_8.sh [graph_abbr]
```

**8. Generate data for Table 9:**
```
bash ./scripts/extract_table_9.sh [graph_abbr]
```
The cache misses are collected by `perf` tool, please make sure the CPU supports related hardware counters.

**9. Generate data for Table 10:**
```
bash ./scripts/extract_table_10.sh [graph_abbr]
```

**10. Generate data for Table 12:**
```
bash ./scripts/extract_table_12.sh [graph_abbr]
```

**11. Generate data for Table 13:**
```
../apps/SSSP_Batch -option ground-truth -batch 2 -max_combination 512 -qf ../query_input/LJ_queries.txt ../ae-data/soc-LiveJournal1.weighted.adj
```
The results table will be printed in the terminal after the program finishes (estimation: around 1 hour).

**12. Data for Table 14:**
Please check any of Glign's output file inside the results folder, there is a line `Profiling cost: ` in each file.
For example, in `results/LJ/BFS_glign_batch_LJ_2.txt` line 7, the profiling cost is 0.34739 seconds.

**13. Data for Table 15:**
```
bash ./scripts/extract_table_15.sh [graph_abbr]
```
| **Graph** | **Command** |
|-------------|------------------------------|
| RD-CA | bash ./scripts/extract_table_15.sh RDCA |
| RD-US | bash ./scripts/extract_table_15.sh RDUS |

**14. Data for Table 16:**
```
bash ./scripts/extract_table_16.sh [graph_abbr]
```
