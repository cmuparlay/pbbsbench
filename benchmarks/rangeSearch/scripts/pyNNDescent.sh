#!/bin/bash
cd ~/pbbsbench/benchmarks/rangeSearch/pyNNDescent

make clean all
P=/ssd1/data/FB_ssnpp
R=/ssd1/results/FB_ssnpp
./range -R 60 -L 1000 -a 10 -d 1.2 -b 2 -q $P/FB_ssnpp_public_queries.u8bin -o $R/1M_60 -c $P/ssnpp-1M -res $R/pynn_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_1000000
./range -R 90 -L 1000 -a 10 -d 1.2 -b 2 -q $P/FB_ssnpp_public_queries.u8bin -o $R/1M_90 -c $P/ssnpp-1M -res $R/pynn_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_1000000

./range -R 60 -L 1000 -a 10 -d 1.2 -b 2 -q $P/FB_ssnpp_public_queries.u8bin -o $R/10M_60 -c $P/ssnpp-10M -res $R/pynn_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_10000000
./range -R 90 -L 1000 -a 10 -d 1.2 -b 2 -q $P/FB_ssnpp_public_queries.u8bin -o $R/10M_90 -c $P/ssnpp-10M -res $R/pynn_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_10000000

./range -R 60 -L 1000 -a 10 -d 1.2 -b 2 -q $P/FB_ssnpp_public_queries.u8bin -o $R/100M_60 -c $P/ssnpp-100M -res $R/pynn_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_100000000
./range -R 90 -L 1000 -a 10 -d 1.2 -b 2 -q $P/FB_ssnpp_public_queries.u8bin -o $R/100M_90 -c $P/ssnpp-100M -res $R/pynn_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_100000000