#!/bin/bash
cd ~/pbbsbench/benchmarks/ANN/pyNNDescent
make clean all

P=/ssd1/data
G=/ssd1/results
#BIGANN: 1M runs
BP=$P/bigann
BG=$G/bigann
./neighbors -R 40 -L 1000 -a 10 -d 1.2 -b 2 -o $BG/1M_pynn_40 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/pynn_res_new.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
./neighbors -R 80 -L 1000 -a 10 -d 1.2 -b 2 -o $BG/1M_pynn_80 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/pynn_res_new.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000

SP=$P/MSSPACEV1B
SG=$G/MSSPACEV1B
#SPACEV:10M and 100M, queries only
./neighbors -R 90 -L 1000 -a 10 -d 1.2 -b 2 -g $SG/10M_pynn_90 -q $SP/query.i8bin -c $SP/msspacev-10M -res $SG/pynn_res_new.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_10000000
./neighbors -R 90 -L 1000 -a 10 -d 1.2 -b 2 -g $SG/100M_pynn_90 -q $SP/query.i8bin -c $SP/msspacev-100M -res $SG/pynn_res_new.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_100000000
