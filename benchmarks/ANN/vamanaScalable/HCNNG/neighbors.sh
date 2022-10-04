#!/bin/bash

make 

# ./neighbors -a 1000 -R 3 -L 10 -k 200 -Q 250 -b 1 -q /ssd1/ANN/sift/bigann_query.bvecs -o outFile -c /ssd1/ANN/sift/idx_1M.ivecs /ssd1/ANN/sift/sift1M.bvecs
# ./neighbors -a 1000 -R 3 -L 30 -k 200 -Q 250 -b 1 -q /ssd1/ANN/sift/bigann_query.bvecs -o outFile -c /ssd1/ANN/sift/idx_1M.ivecs /ssd1/ANN/sift/sift1M.bvecs
# ./neighbors -a 1000 -R 3 -L 10 -k 200 -Q 250 -b 1 -q /ssd1/ANN/sift/bigann_query.bvecs -o outFile -c /ssd1/ANN/sift/idx_10M.ivecs /ssd1/ANN/sift/sift10M.bvecs
./neighbors -a 1000 -R 3 -L 30 -k 200 -Q 250 -b 1 -q /ssd1/ANN/sift/bigann_query.bvecs -o outFile -c /ssd1/ANN/sift/idx_10M.ivecs /ssd1/ANN/sift/sift10M.bvecs
# ./neighbors -a 1000 -R 3 -L 10 -k 200 -Q 250 -b 1 -q /ssd1/ANN/sift/bigann_query.bvecs -o outFile -c /ssd1/ANN/sift/idx_1000M.ivecs /ssd1/ANN/sift/bigann_base.bvecs
# ./neighbors -a 1000 -R 3 -L 30 -k 200 -Q 250 -b 1 -q /ssd1/ANN/sift/bigann_query.bvecs -o outFile -c /ssd1/ANN/sift/idx_1000M.ivecs /ssd1/ANN/sift/bigann_base.bvecs