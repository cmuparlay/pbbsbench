#!/bin/bash

make
./range -a 1000 -R 3 -L 10 -k 200 -Q 250 -b 1 -q /ssd1/bigann/FB_ssnpp/FB_ssnpp_public_queries.u8bin -o outFile -c /ssd1/bigann/FB_ssnpp/ssnpp-10M /ssd1/bigann/FB_ssnpp/FB_ssnpp_database.u8bin.crop_nb_10000000