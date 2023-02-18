#!/bin/bash
cd ~/pbbsbench/benchmarks/rangeSearch/HCNNG

make clean all
P=/ssd1/data/FB_ssnpp
R=/ssd1/results/FB_ssnpp
./range -a 1000 -R 3 -L 30 -b 1 -q $P/FB_ssnpp_public_queries.u8bin -o $R/1M_3_30 -c $P/ssnpp-1M -res $R/hcnng_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_1000000
./range -a 1000 -R 3 -L 50 -b 1 -q $P/FB_ssnpp_public_queries.u8bin -o $R/1M_3_50 -c $P/ssnpp-1M -res $R/hcnng_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_1000000

./range -a 1000 -R 3 -L 30 -b 1 -q $P/FB_ssnpp_public_queries.u8bin -o $R/10M_3_30 -c $P/ssnpp-10M -res $R/hcnng_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_10000000
./range -a 1000 -R 3 -L 50 -b 1 -q $P/FB_ssnpp_public_queries.u8bin -o $R/10M_3_50 -c $P/ssnpp-10M -res $R/hcnng_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_10000000

./range -a 1000 -R 3 -L 30 -b 1 -q $P/FB_ssnpp_public_queries.u8bin -o $R/100M_3_30 -c $P/ssnpp-100M -res $R/hcnng_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_100000000
./range -a 1000 -R 3 -L 50 -b 1 -q $P/FB_ssnpp_public_queries.u8bin -o $R/100M_3_50 -c $P/ssnpp-100M -res $R/hcnng_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_100000000

./range -a 1000 -R 3 -L 50 -b 1 -q $P/FB_ssnpp_public_queries.u8bin -o $R/1B_3_50 -c $P/FB_ssnpp_public_queries_1B_GT.rangeres -res $R/hcnng_res.csv $P/FB_ssnpp_database.u8bin