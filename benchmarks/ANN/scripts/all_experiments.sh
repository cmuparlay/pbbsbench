#!/bin/bash
#missing queries
mkdir /ssd1/results/piped_output
cd ~/pbbsbench/benchmarks/ANN/scripts
bash vamana.sh
bash hcnng.sh
bash pynndescent.sh
cd ~/pbbsbench/benchmarks/rangeSearch/scripts
bash search_all.sh
#FAISS with better strings for all datasets: est 16 hours
cd ~/pbbsbench/benchmarks/ANN/scripts/FAISS
bash run2.sh



