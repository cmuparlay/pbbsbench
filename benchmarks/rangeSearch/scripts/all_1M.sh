#!/bin/bash
P=/ssd1/data/FB_ssnpp
R=/ssd1/results/FB_ssnpp

cd ~/pbbsbench/benchmarks/rangeSearch/vamana
make clean all
./range -a 1.2 -R 128 -L 256 -q $P/FB_ssnpp_public_queries.u8bin -o $R/1M_128_256 -c $P/ssnpp-1M -res $R/1M_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_1000000

cd ~/pbbsbench/benchmarks/rangeSearch/HCNNG
make clean all
./range -a 1000 -R 3 -L 30 -b 1 -q $P/FB_ssnpp_public_queries.u8bin -o $R/1M_3_30 -c $P/ssnpp-1M -res $R/hcnng_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_1000000

cd ~/pbbsbench/benchmarks/rangeSearch/pyNNDescent
make clean all
./range -R 60 -L 1000 -a 10 -d 1.2 -b 2 -q $P/FB_ssnpp_public_queries.u8bin -o $R/1M_60 -c $P/ssnpp-1M -res $R/pynn_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_1000000