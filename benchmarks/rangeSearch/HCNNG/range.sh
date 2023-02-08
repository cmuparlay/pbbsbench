#!/bin/bash

make
./range -a 1000 -R 3 -L 30 -b 1 -q /ssd1/bigann/FB_ssnpp/FB_ssnpp_public_queries.u8bin -res hcnng_res.csv -c /ssd1/bigann/FB_ssnpp/gt_1M /ssd1/bigann/FB_ssnpp/FB_ssnpp_database.u8bin.crop_nb_1000000
./range -a 1000 -R 3 -L 50 -b 1 -q /ssd1/bigann/FB_ssnpp/FB_ssnpp_public_queries.u8bin -res hcnng_res.csv -c /ssd1/bigann/FB_ssnpp/gt_1M /ssd1/bigann/FB_ssnpp/FB_ssnpp_database.u8bin.crop_nb_1000000

./range -a 1000 -R 3 -L 30 -b 1 -q /ssd1/bigann/FB_ssnpp/FB_ssnpp_public_queries.u8bin -res hcnng_res.csv -c /ssd1/bigann/FB_ssnpp/ssnpp-10M /ssd1/bigann/FB_ssnpp/FB_ssnpp_database.u8bin.crop_nb_10000000
./range -a 1000 -R 3 -L 50 -b 1 -q /ssd1/bigann/FB_ssnpp/FB_ssnpp_public_queries.u8bin -res hcnng_res.csv -c /ssd1/bigann/FB_ssnpp/ssnpp-10M /ssd1/bigann/FB_ssnpp/FB_ssnpp_database.u8bin.crop_nb_10000000

./range -a 1000 -R 3 -L 30 -b 1 -q /ssd1/bigann/FB_ssnpp/FB_ssnpp_public_queries.u8bin -res hcnng_res.csv -c /ssd1/bigann/FB_ssnpp/ssnpp-100M /ssd1/bigann/FB_ssnpp/FB_ssnpp_database.u8bin.crop_nb_100000000
./range -a 1000 -R 3 -L 50 -b 1 -q /ssd1/bigann/FB_ssnpp/FB_ssnpp_public_queries.u8bin -res hcnng_res.csv -c /ssd1/bigann/FB_ssnpp/ssnpp-100M /ssd1/bigann/FB_ssnpp/FB_ssnpp_database.u8bin.crop_nb_100000000