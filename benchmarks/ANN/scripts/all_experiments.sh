#!/bin/bash
#TODO separate based on machine?
#minimal preparation: est 10 mins
cd ~/pbbsbench/benchmarks/ANN/scripts
bash prepare_datasets.sh
#pynndescent on all datasets up to 100M: est 12 hours
cd ~/pbbsbench/benchmarks/ANN/scripts
bash pynndescent.sh
#FAISS with better strings for all datasets: est 16 hours
cd ~/big-ann-benchmarks
rm results
mkdir /ssd1/results/FAISSresults2
ln -s /ssd1/results/FAISSresults2 results
cd ~/pbbsbench/benchmarks/ANN/scripts/FAISS
bash run2.sh
#HCNNG: mostly queries, one billion scale build: est 24 h
cd ~/pbbsbench/benchmarks/ANN/scripts
bash hcnng.sh
#range search: 1 billion scale build with HCNNG, plus searches: est 18 hours
#TODO check whether HCNNG build went through
cd ~/pbbsbench/benchmarks/rangeSearch/scripts
bash search_all.sh
#vamana: 2 billion scale builds, with low parameters: est 
cd ~/pbbsbench/benchmarks/ANN/scripts
bash vamana.sh


