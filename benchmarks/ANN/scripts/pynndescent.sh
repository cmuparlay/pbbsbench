#!/bin/bash
cd ~/pbbsbench/benchmarks/ANN/pyNNDescent
make

P=/ssd1/data
G=/ssd1/results
#BIGANN: two settings
BP=$P/bigann
BG=$G/bigann
./neighbors -R 40 -L 1000 -a 10 -d 1.2 -b 2 -o $BG/1M_pynn_40 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/pynn_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
./neighbors -R 40 -L 1000 -a 10 -d 1.2 -b 2 -o $BG/10M_pynn_40 -q $BP/query.public.10K.u8bin -c $BP/bigann-10M -res $BG/pynn_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_10000000
./neighbors -R 40 -L 1000 -a 10 -d 1.2 -b 2 -o $BG/100M_pynn_40 -q $BP/query.public.10K.u8bin -c $BP/bigann-100M -res $BG/pynn_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_100000000

./neighbors -R 80 -L 1000 -a 10 -d 1.2 -b 2 -o $BG/1M_pynn_80 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/pynn_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
./neighbors -R 80 -L 1000 -a 10 -d 1.2 -b 2 -o $BG/10M_pynn_80 -q $BP/query.public.10K.u8bin -c $BP/bigann-10M -res $BG/pynn_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_10000000
./neighbors -R 80 -L 1000 -a 10 -d 1.2 -b 2 -o $BG/100M_pynn_80 -q $BP/query.public.10K.u8bin -c $BP/bigann-100M -res $BG/pynn_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_100000000
#MSSPACEV: two settings
SP=$P/MSSPACEV1B
SG=$G/MSSPACEV1B
./neighbors -R 60 -L 1000 -a 10 -d 1.2 -b 2 -o $SG/1M_pynn_60 -q $SP/query.i8bin -c $SP/msspacev-1M -res $SG/pynn_res.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_1000000
./neighbors -R 60 -L 1000 -a 10 -d 1.2 -b 2 -o $SG/10M_pynn_60 -q $SP/query.i8bin -c $SP/msspacev-10M -res $SG/pynn_res.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_10000000
./neighbors -R 60 -L 1000 -a 10 -d 1.2 -b 2 -o $SG/100M_pynn_60 -q $SP/query.i8bin -c $SP/msspacev-100M -res $SG/pynn_res.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_100000000

./neighbors -R 90 -L 1000 -a 10 -d 1.2 -b 2 -o $SG/1M_pynn_90 -q $SP/query.i8bin -c $SP/msspacev-1M -res $SG/pynn_res.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_1000000
./neighbors -R 90 -L 1000 -a 10 -d 1.2 -b 2 -o $SG/10M_pynn_90 -q $SP/query.i8bin -c $SP/msspacev-10M -res $SG/pynn_res.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_10000000
./neighbors -R 90 -L 1000 -a 10 -d 1.2 -b 2 -o $SG/100M_pynn_90 -q $SP/query.i8bin -c $SP/msspacev-100M -res $SG/pynn_res.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_100000000
# #TEXT2IMAGE: two settings
TP=$P/text2image1B
TG=$G/text2image1B
./neighbors -R 60 -L 1000 -a 10 -d 1.2 -b 2 -o $TG/1M_pynn_90 -q $TP/query.public.100K.fbin -c $TP/text2image-1M -res $TG/pynn_res.csv -f bin -t float $TP/base.1B.fbin.crop_nb_1000000
./neighbors -R 60 -L 1000 -a 10 -d 1.2 -b 2 -o $TG/10M_pynn_90 -q $TP/query.public.100K.fbin -c $TP/text2image-10M -res $TG/pynn_res.csv -f bin -t float $TP/base.1B.fbin.crop_nb_10000000
./neighbors -R 60 -L 1000 -a 10 -d 1.2 -b 2 -o $TG/100M_pynn_90 -q $TP/query.public.100K.fbin -c $TP/text2image-100M -res $TG/pynn_res.csv -f bin -t float $TP/base.1B.fbin.crop_nb_100000000

./neighbors -R 90 -L 1000 -a 10 -d 1.2 -b 2 -o $TG/1M_pynn_90 -q $TP/query.public.100K.fbin -c $TP/text2image-1M -res $TG/pynn_res.csv -f bin -t float $TP/base.1B.fbin.crop_nb_1000000
./neighbors -R 90 -L 1000 -a 10 -d 1.2 -b 2 -o $TG/10M_pynn_90 -q $TP/query.public.100K.fbin -c $TP/text2image-10M -res $TG/pynn_res.csv -f bin -t float $TP/base.1B.fbin.crop_nb_10000000
./neighbors -R 90 -L 1000 -a 10 -d 1.2 -b 2 -o $TG/100M_pynn_90 -q $TP/query.public.100K.fbin -c $TP/text2image-100M -res $TG/pynn_res.csv -f bin -t float $TP/base.1B.fbin.crop_nb_100000000

