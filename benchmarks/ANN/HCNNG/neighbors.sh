#!/bin/bash

make 
# ./neighbors -a 1000 -R 3 -L 20 -k 200 -Q 250 -b 1 -o /ssd1/ANN/sift/1M_3_20 /ssd1/ANN/sift/sift1M.bvecs
# ./neighbors -a 1000 -R 3 -L 20 -k 200 -Q 250 -b 1 -q /ssd1/ANN/sift/bigann_query.bvecs -g /ssd1/ANN/sift/1M_3_20 -c /ssd1/ANN/sift/idx_1M.ivecs /ssd1/ANN/sift/sift1M.bvecs
./neighbors -a 1000 -R 3 -L 20 -k 200 -Q 250 -b 1 -q /ssd1/ANN/sift/bigann_query.bvecs -c /ssd1/ANN/sift/idx_1M.ivecs -res hcnng_res.csv -f vec -t uint8 /ssd1/ANN/sift/sift1M.bvecs

# make clean all
# echo ""
# echo "10 trees"
# ./neighbors -a 1000 -R 3 -L 10 -k 200 -Q 250 -b 1 -q /ssd1/ANN/yandex/yandex_query.fvecs -o outFile -c /ssd1/ANN/yandex/yandex_groundtruth.ivecs /ssd1/ANN/yandex/yandex_1M.fvecs
# echo ""
# echo "20 trees"
# ./neighbors -a 1000 -R 3 -L 20 -k 200 -Q 250 -b 1 -q /ssd1/ANN/yandex/yandex_query.fvecs -o outFile -c /ssd1/ANN/yandex/yandex_groundtruth.ivecs /ssd1/ANN/yandex/yandex_1M.fvecs
# echo ""
# echo "30 trees"
# ./neighbors -a 1000 -R 3 -L 30 -k 200 -Q 250 -b 1 -q /ssd1/ANN/yandex/yandex_query.fvecs -o outFile -c /ssd1/ANN/yandex/yandex_groundtruth.ivecs /ssd1/ANN/yandex/yandex_1M.fvecs
# echo ""
# echo "30 trees, 1B points"
# ./neighbors -a 1000 -R 3 -L 30 -k 200 -Q 250 -b 1 -q /ssd1/ANN/sift/bigann_query.bvecs -o outFile -c /ssd1/ANN/sift/idx_1000M.ivecs /ssd1/ANN/sift/bigann_base.bvecs
# ./neighbors -a 1000 -R 3 -L 10 -k 200 -Q 250 -b 1 -q /ssd1/ANN/sift/bigann_query.bvecs -o outFile -c /ssd1/ANN/sift/idx_1000M.ivecs /ssd1/ANN/sift/bigann_base.bvecs
# ./neighbors -a 1000 -R 3 -L 30 -k 200 -Q 250 -b 1 -q /ssd1/ANN/sift/bigann_query.bvecs -o outFile -c /ssd1/ANN/sift/idx_1000M.ivecs /ssd1/ANN/sift/bigann_base.bvecs