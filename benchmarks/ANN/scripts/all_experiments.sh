#!/bin/bash
cd ~/pbbsbench/benchmarks/ANN/scripts
#download datasets
# bash download_datasets.sh
#prepare groundtruth for 1M slices
# bash prepare_datasets.sh
#nearest neighbor search experiments
nohup bash hcnng.sh >> run_log.out
nohup bash vamana.sh >> run_log.out
# bash pynndescent.sh
#range search experiments
# cd ~/pbbsbench/benchmarks/rangeSearch/scripts
# bash all_experiments.sh
#all experiments with FAISS
# cd ~/pbbsbench/benchmarks/ANN/scripts
# bash run.sh



