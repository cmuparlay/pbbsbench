#/bin/bash
cd ~/hcnng
D=sift/sift_base.fvecs
P=~/pbbsbench/benchmarks/ANN/scripts/parallel_tests
g++ hcnng.cpp -o hcnng -std=c++11 -fopenmp -O3
for i in 1 2 8 12 24 30; do
    export OMP_NUM_THREADS=$i
    ./hcnng $D 1000 30 hcnng_sift.ivecs >> $P/hcnng_data.txt
done
