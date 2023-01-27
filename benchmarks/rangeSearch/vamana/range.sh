#!/bin/bash

make
# echo "1M, R=64"
./range -a 1.2 -R 64 -L 128 -k 10 -Q 20 -q /ssd1/bigann/FB_ssnpp/FB_ssnpp_public_queries.u8bin -o /ssd1/bigann/FB_ssnpp/1M_64_128 -c /ssd1/bigann/FB_ssnpp/gt_1M /ssd1/bigann/FB_ssnpp/FB_ssnpp_database.u8bin.crop_nb_1000000
# echo "1M, R=128"
# ./range -a 1.2 -R 128 -L 256 -k 10 -Q 20 -q /ssd1/bigann/FB_ssnpp/FB_ssnpp_public_queries.u8bin -o /ssd1/bigann/FB_ssnpp/1M_128_256 -c /ssd1/bigann/FB_ssnpp/gt_1M /ssd1/bigann/FB_ssnpp/FB_ssnpp_database.u8bin.crop_nb_1000000
# echo "10M, R=64"
# ./range -a 1.2 -R 64 -L 128 -k 10 -Q 20 -q /ssd1/bigann/FB_ssnpp/FB_ssnpp_public_queries.u8bin -g /ssd1/bigann/FB_ssnpp/10M_64_128 -c /ssd1/bigann/FB_ssnpp/ssnpp-10M /ssd1/bigann/FB_ssnpp/FB_ssnpp_database.u8bin.crop_nb_10000000
# echo "10M, R=128"
# ./range -a 1.2 -R 128 -L 256 -k 10 -Q 20 -q /ssd1/bigann/FB_ssnpp/FB_ssnpp_public_queries.u8bin -o /ssd1/bigann/FB_ssnpp/10M_128_256 -c /ssd1/bigann/FB_ssnpp/ssnpp-10M /ssd1/bigann/FB_ssnpp/FB_ssnpp_database.u8bin.crop_nb_10000000
# ./range -a 1.2 -R 64 -L 128 -k 50 -Q 100 -q /ssd1/bigann/FB_ssnpp/FB_ssnpp_public_queries.u8bin -g /ssd1/bigann/FB_ssnpp/1M_64_128 -c /ssd1/bigann/FB_ssnpp/gt_1M /ssd1/bigann/FB_ssnpp/FB_ssnpp_database.u8bin.crop_nb_1000000