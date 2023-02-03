#!/bin/bash
cd ~/pbbsbench/benchmarks/ANN/vamana
make clean all

P=/ssd1/data
G=/ssd1/results
# BIGANN
BP=$P/bigann
BG=$G/bigann
./neighbors -R 64 -L 128 -o $BG/1M_vamana_64_128 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/vamana_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
#MSSPACEV
SP=$P/MSSPACEV1B
SG=$G/MSSPACEV1B
./neighbors -R 64 -L 128 -o $SG/1M_vamana_64_128 -q $SP/query.i8bin -c $SP/msspacev-1M -res $SG/vamana_res.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_1000000
#TEXT2IMAGE
TP=$P/text2image1B
TG=$G/text2image1B
./neighbors -R 64 -L 128 -o $TG/1M_vamana_64_128 -q $TP/query.public.100K.fbin -c $TP/text2image-1M -res $TG/vamana_res.csv -f bin -t float $TP/base.1B.fbin.crop_nb_1000000

#start HCNNG
cd ~/pbbsbench/benchmarks/ANN/HCNNG
make
#BIGANN
./neighbors -a 1000 -R 3 -L 20 -b 1 -o $BG/1M_HCNNG_20 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/HCNNG_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
#MSSPACEV:
./neighbors -a 1000 -R 3 -L 30 -b 1 -o $SG/1M_HCNNG_30 -q $SP/query.i8bin -c $SP/msspacev-1M -res $SG/HCNNG_res.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_1000000
# #TEXT2IMAGE:
./neighbors -a 1000 -R 3 -L 30 -b 1 -o $TG/1M_HCNNG_30 -q $TP/query.public.100K.fbin -c $TP/text2image-1M -res $TG/HCNNG_res.csv -f bin -t float $TP/base.1B.fbin.crop_nb_1000000

#start PyNNDescent
cd ~/pbbsbench/benchmarks/ANN/pyNNDescent
make
#BIGANN
./neighbors -R 40 -L 1000 -a 10 -d 1.2 -b 2 -o $BG/1M_pynn_40 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/pynn_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
#MSSPACEV
./neighbors -R 60 -L 1000 -a 10 -d 1.2 -b 2 -o $SG/1M_pynn_60 -q $SP/query.i8bin -c $SP/msspacev-1M -res $SG/pynn_res.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_1000000
# #TEXT2IMAGE
./neighbors -R 60 -L 1000 -a 10 -d 1.2 -b 2 -o $TG/1M_pynn_90 -q $TP/query.public.100K.fbin -c $TP/text2image-1M -res $TG/pynn_res.csv -f bin -t float $TP/base.1B.fbin.crop_nb_1000000

