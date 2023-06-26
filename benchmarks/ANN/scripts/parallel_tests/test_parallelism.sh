#!/bin/bash

BP=/ssd1/data/bigann
G=/ssd1/results/bigann
D=~/pbbsbench/benchmarks/ANN/scripts/parallel_tests

# echo "test run"
# cd ~/pbbsbench/benchmarks/ANN/vamana
# make clean all 
# ./neighbors -R 64 -L 128 -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
# cd ~/pbbsbench/benchmarks/ANN/HCNNG
# make clean all 
# ./neighbors -a 1000 -R 3 -L 30 -b 1 -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
# cd ~/pbbsbench/benchmarks/ANN/pyNNDescent
# make clean all 
# ./neighbors -R 60 -L 100 -a 10 -d 1.2 -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000



cd ~/pbbsbench/benchmarks/ANN/vamana
make 
for i in 1 2 8 12 24 48 96 144 192; do
    PARLAY_NUM_THREADS=$i nohup ./neighbors -R 64 -L 128 -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000 >> $D/test_parallel.out
done

cd ~/pbbsbench/benchmarks/ANN/HCNNG
make 
for i in 1 2 8 24 48 96 144 192; do
    PARLAY_NUM_THREADS=$i nohup ./neighbors -a 1000 -R 3 -L 30 -b 1 -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000 >> $D/test_parallel.out
done

cd ~/pbbsbench/benchmarks/ANN/pyNNDescent
make 
for i in 1 2 8 24 48 96 144 192; do
    PARLAY_NUM_THREADS=$i nohup ./neighbors -R 60 -L 100 -a 10 -d 1.2 -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000 >> $D/test_parallel.out
done
