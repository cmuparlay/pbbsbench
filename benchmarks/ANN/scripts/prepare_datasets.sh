#calculate groundtruth files for all 1M slices
cd ~/pbbsbench/benchmarks/ANN/vamana
make compute_groundtruth
P=/ssd1/data
TP=$P/text2image1B
./compute_groundtruth $TP/base.1B.fbin.crop_nb_1000000 $TP/query.public.100K.fbin "bin" "float" 100 1 $TP/text2image-1M
