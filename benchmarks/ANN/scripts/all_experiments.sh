#!/bin/bash
cd ~/pbbsbench/benchmarks/ANN/scripts
#download datasets
bash download_datasets.sh
#prepare groundtruth for 1M slices
bash prepare_datasets.sh
#nearest neighbor search experiments
bash vamana.sh
bash hcnng.sh
bash pynndescent.sh
#range search experiments
cd ~/pbbsbench/benchmarks/rangeSearch/scripts
bash all_experiments.sh
#all experiments with FAISS
cd ~/pbbsbench/benchmarks/ANN/scripts
bash run.sh



