#!/bin/bash
cd ~/pbbsbench/benchmarks/ANN/HCNNG
make clean all

P=/ssd1/data
G=/ssd1/results
BP=$P/bigann
BG=$G/bigann

#test with varying num_threads
PARLAY_NUM_THREADS=1 ./neighbors -a 1000 -R 3 -L 10 -b 1 -o $BG/1M_HCNNG_30 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/test_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
PARLAY_NUM_THREADS=2 ./neighbors -a 1000 -R 3 -L 10 -b 1 -o $BG/1M_HCNNG_30 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/test_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
PARLAY_NUM_THREADS=48 ./neighbors -a 1000 -R 3 -L 10 -b 1 -o $BG/1M_HCNNG_30 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/test_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
PARLAY_NUM_THREADS=96 ./neighbors -a 1000 -R 3 -L 10 -b 1 -o $BG/1M_HCNNG_30 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/test_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
./neighbors -a 1000 -R 3 -L 30 -b 1 -o $BG/1M_HCNNG_30 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/test_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000

#test with varying num_threads and numa
PARLAY_NUM_THREADS=1 numactl -i all ./neighbors -a 1000 -R 3 -L 10 -b 1 -o $BG/1M_HCNNG_30 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/test_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
PARLAY_NUM_THREADS=2 numactl -i all ./neighbors -a 1000 -R 3 -L 10 -b 1 -o $BG/1M_HCNNG_30 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/test_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
PARLAY_NUM_THREADS=48 numactl -i all  ./neighbors -a 1000 -R 3 -L 10 -b 1 -o $BG/1M_HCNNG_30 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/test_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
PARLAY_NUM_THREADS=96 numactl -i all  ./neighbors -a 1000 -R 3 -L 10 -b 1 -o $BG/1M_HCNNG_30 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/test_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
numactl -i all ./neighbors -a 1000 -R 3 -L 10 -b 1 -o $BG/1M_HCNNG_30 -q $BP/query.public.10K.u8bin -c $BP/bigann-1M -res $BG/test_res.csv -f bin -t uint8 $BP/base.1B.u8bin.crop_nb_1000000
