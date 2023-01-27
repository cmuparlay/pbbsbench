#!/bin/bash

make 
./range -R 40 -L 1000 -a 10 -d 1.2 -Q 250 -k 200 -b 2 -q /ssd1/bigann/FB_ssnpp/FB_ssnpp_public_queries.u8bin -o outFile -c /ssd1/bigann/FB_ssnpp/ssnpp-10M /ssd1/bigann/FB_ssnpp/FB_ssnpp_database.u8bin.crop_nb_10000000