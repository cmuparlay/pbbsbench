#!/bin/bash
cd ~/pbbsbench/benchmarks/ANN/vamana
make clean all

P=/ssd1/data
G=/ssd1/results
# BIGANN
BP=$P/bigann
BG=$G/bigann
./neighbors -R 64 -L 128 -o $BG/1M_vamana_64_128 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $G/million_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
#MSSPACEV
SP=$P/MSSPACEV1B
SG=$G/MSSPACEV1B
./neighbors -R 64 -L 128 -o $SG/1M_vamana_64_128 -q $SP/query.i8bin -c $SP/msspacev-1M -res $G/million_res.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_1000000
#TEXT2IMAGE
TP=$P/text2image1B
TG=$G/text2image1B
./neighbors -a .9 -R 64 -L 128 -o $TG/1M_vamana_64_128 -q $TP/query.public.100K.fbin -c $TP/text2image-1M -res $G/million_res.csv -f bin -t float -D 1 $TP/base.1B.fbin.crop_nb_1000000

#start HCNNG
cd ~/pbbsbench/benchmarks/ANN/HCNNG
make clean all
#BIGANN
./neighbors -a 1000 -R 3 -L 30 -b 1 -o $BG/1M_HCNNG_30 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $O/bigann_HCNNG.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
#MSSPACEV:
./neighbors -a 1000 -R 3 -L 50 -b 1 -o $SG/1M_HCNNG_50 -q $SP/query.i8bin -c $SP/msspacev-1M -res $O/spacev_HCNNG.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_1000000
# #TEXT2IMAGE:
./neighbors -a 1000 -R 3 -L 30 -b 1 -o $TG/1M_HCNNG_30 -q $TP/query.public.100K.fbin -c $TP/text2image-1M -res $TG/t2i_HCNNG.csv -f bin -t float -D 1 $TP/base.1B.fbin.crop_nb_1000000

#start PyNNDescent
cd ~/pbbsbench/benchmarks/ANN/pyNNDescent
make clean all
#BIGANN
./neighbors -R 40 -L 100 -a 10 -d 1.2 -o $BG/1M_pynn_40 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $G/bigann_pynn.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
#MSSPACEV
./neighbors -R 60 -L 100 -a 10 -d 1.2 -o $SG/1M_pynn_60 -q $SP/query.i8bin -c $SP/msspacev-1M -res $G/spacev_pynn.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_1000000
# #TEXT2IMAGE
./neighbors -R 60 -L 100 -a 10 -d .9 -o $TG/1M_pynn_60 -q $TP/query.public.100K.fbin -c $TP/text2image-1M -res $G/t2i_pynn.csv -f bin -t float -D 1 $TP/base.1B.fbin.crop_nb_1000000

