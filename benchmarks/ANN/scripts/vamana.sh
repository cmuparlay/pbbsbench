#!/bin/bash
cd ~/pbbsbench/benchmarks/ANN/vamana
make clean all

P=/ssd1/data
G=/ssd1/results
BP=$P/bigann
BG=$G/bigann
# ./neighbors -R 64 -L 128 -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
# ./neighbors -R 64 -L 128 -o $BG/1M_vamana_64_128 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $G/bigann_vamana.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
# ./neighbors -R 64 -L 128 -o $BG/10M_vamana_64_128 -q $BP/query.public.10K.u8bin -c $BP/bigann-10M -res $G/bigann_vamana.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_10000000
# ./neighbors -R 64 -L 128 -o $BG/100M_vamana_64_128 -q $BP/query.public.10K.u8bin -c $BP/bigann-100M -res $G/bigann_vamana.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_100000000
# ./neighbors -R 64 -L 128 -o $BG/1B_vamana_64_128 -q $BP/query.public.10K.u8bin -c $BP/bigann-1B -res $G/bigann_vamana.csv -f bin -t uint8 $BP/base.1B.u8bin


SP=$P/MSSPACEV1B
SG=$G/MSSPACEV1B
# ./neighbors -R 64 -L 128 -o $SG/1M_vamana_64_128 -q $SP/query.i8bin -c $SP/msspacev-1M -res $G/spacev_vamana.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_1000000
# ./neighbors -R 64 -L 128 -o $SG/10M_vamana_64_128 -q $SP/query.i8bin -c $SP/msspacev-10M -res $G/spacev_vamana.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_10000000
# ./neighbors -R 64 -L 128 -o $SG/100M_vamana_64_128 -q $SP/query.i8bin -c $SP/msspacev-100M -res $G/spacev_vamana.csv -f bin -t int8 $SP/spacev1b_base.i8bin.crop_nb_100000000
# ./neighbors -R 64 -L 128 -o $SG/1B_vamana_64_128 -q $SP/query.i8bin -c $SP/msspacev-1B -res $G/spacev_vamana.csv -f bin -t int8 $SP/spacev1b_base.i8bin


TP=$P/text2image1B
TG=$G/text2image1B
# ./neighbors -a .9 -R 64 -L 128 -o $TG/1M_vamana_64_128 -q $TP/query.public.100K.fbin -c $TP/text2image-1M -res $G/spacev_vamana.csv -f bin -t float -D 1 $TP/base.1B.fbin.crop_nb_1000000
# ./neighbors -a .9 -R 64 -L 128 -o $TG/10M_vamana_64_128 -q $TP/query.public.100K.fbin -c $TP/text2image-10M -res $G/spacev_vamana.csv -f bin -t float -D 1 $TP/base.1B.fbin.crop_nb_10000000
# ./neighbors -a .98 -R 64 -L 128 -o $TG/100M_vamana_64_128 -q $TP/query.public.100K.fbin -c $TP/text2image-100M -res $G/t2i_vamana_resubmit.csv -f bin -t float -D 1 $TP/base.1B.fbin.crop_nb_100000000
./neighbors -a .98 -R 64 -L 128 -o $TG/1B_vamana_64_128 -q $TP/query.public.100K.fbin -c $TP/text2image-1B -res $G/t2i_vamana_resubmit.csv -f bin -t float -D 1 $TP/base.1B.fbin

