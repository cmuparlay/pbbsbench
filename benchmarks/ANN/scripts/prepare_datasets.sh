#calculate groundtruth files for all 1M slices
cd ~/pbbsbench/benchmarks/ANN/vamana
make compute_groundtruth
make crop_sift
P=/ssd1/data
BP=$P/bigann
./crop_sift $BP/base.1B.u8bin.crop_nb_10000000 $BP/base.1B.u8bin.crop_nb_1000000
./compute_groundtruth $BP/base.1B.u8bin.crop_nb_1000000 $BP/query.public.10K.u8bin "bin" "uint8" 100 0 $BP/bigann-1M
SP=$P/MSSPACEV1B
./compute_groundtruth $SP/spacev1b_base.i8bin.crop_nb_1000000 $SP/query.i8bin "bin" "int8" 100 0 $SP/msspacev-1M
TP=$P/text2image1B
./compute_groundtruth $TP/base.1B.fbin.crop_nb_1000000 $TP/query.public.100K.fbin "bin" "float" 100 1 $TP/text2image-1M
cd ../../rangeSearch/vamana
SP=$P/FB_ssnpp
make compute_range_groundtruth
./compute_range_groundtruth $SP/FB_ssnpp_database.u8bin.crop_nb_1000000 $SP/FB_ssnpp_public_queries.u8bin 96237 $SP/ssnpp-1M
