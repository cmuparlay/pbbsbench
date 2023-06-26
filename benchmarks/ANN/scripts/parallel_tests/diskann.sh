#/bin/bash
P=/ssd1/data/bigann
D=~/pbbsbench/benchmarks/ANN/scripts/parallel_tests
cd ~/DiskANN/build
for i in 1 2 8 12 24 48 96 144 192; do
    ./tests/build_memory_index  --data_type uint8 --dist_fn l2 --data_path $P/base.1B.u8bin.crop_nb_1000000 --index_path_prefix $P/test_index -R 50 -L 200 --alpha 1.2 -T $i >> $D/diskann_data.txt
done
