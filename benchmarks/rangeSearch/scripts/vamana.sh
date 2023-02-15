#!/bin/bash
cd ~/pbbsbench/benchmarks/rangeSearch/vamana

make clean all
P=/ssd1/data/FB_ssnpp
R=/ssd1/results/FB_ssnpp
./range -a 1.2 -R 128 -L 256 -q $P/FB_ssnpp_public_queries.u8bin -o $R/1M_128_256 -c $P/ssnpp-1M -res $R/vamana_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_1000000
./range -a 1.2 -R 150 -L 400 -q $P/FB_ssnpp_public_queries.u8bin -o $R/1M_150_400 -c $P/ssnpp-1M -res $R/vamana_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_1000000

./range -a 1.2 -R 128 -L 256 -q $P/FB_ssnpp_public_queries.u8bin -o $R/10M_128_256 -c $P/ssnpp-10M -res $R/vamana_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_10000000
./range -a 1.2 -R 150 -L 400 -q $P/FB_ssnpp_public_queries.u8bin -o $R/10M_150_400 -c $P/ssnpp-10M -res $R/vamana_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_10000000

./range -a 1.2 -R 128 -L 256 -q $P/FB_ssnpp_public_queries.u8bin -o $R/100M_128_256 -c $P/ssnpp-100M -res $R/vamana_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_100000000
./range -a 1.2 -R 150 -L 400 -q $P/FB_ssnpp_public_queries.u8bin -o $R/100M_150_400 -c $P/ssnpp-100M -res $R/vamana_res.csv $P/FB_ssnpp_database.u8bin.crop_nb_100000000

./range -a 1.2 -R 150 -L 400 -q $P/FB_ssnpp_public_queries.u8bin -o $R/1B_150_400 -c $P/FB_ssnpp_public_queries_1B_GT.rangeres -res $R/vamana_res.csv $P/FB_ssnpp_database.u8bin