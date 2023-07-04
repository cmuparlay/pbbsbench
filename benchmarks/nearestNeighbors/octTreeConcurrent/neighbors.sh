#!/bin/bash

# in checkNeighbors, r = number of points to calculate exhaustive nbh for

make 

echo "geolife"
numactl -i all ./neighbors -d 3 -k 1 -o oFile /ssd0/geolife/geolife_10M_filtered
cd ../bench
./neighborsCheck -d 3 -k 1 -r 100 /ssd0/geolife/geolife_10M_filtered ../octTreeConcurrent/oFile
cd ../octTree

# echo "2DinCube"
# PARLAY_NUM_THREADS=140 numactl -i all ./neighbors -d 2 -k 1 -o oFile ../geometryData/data/2DinCube_10M
# cd ../bench
# ./neighborsCheck -d 2 -k 1 -r 100 ../geometryData/data/2DinCube_10M ../octTreeConcurrent/oFile
# cd ../octTreeConcurrent
# echo "2DKuzmin"
# numactl -i all ./neighbors -d 2 -o oFile ../geometryData/data/2Dkuzmin_10M
# cd ../bench
# ./neighborsCheck -d 2 -k 1 -r 100 ../geometryData/data/2Dkuzmin_10M ../octTreeConcurrent/oFile
# cd ../octTreeConcurrent
# echo "3DonSphere"
# numactl -i all ./neighbors -d 3 -o oFile ../geometryData/data/3DonSphere_10M
# cd ../bench
# ./neighborsCheck -d 3 -k 1 -r 100 ../geometryData/data/3DonSphere_10M ../octTreeConcurrent/oFile
# cd ../octTreeConcurrent
# echo "3DinCube"
# numactl -i all ./neighbors -d 3 -o oFile ../geometryData/data/3DinCube_10M
# cd ../bench
# ./neighborsCheck -d 3 -k 1 -r 100 ../geometryData/data/3DinCube_10M ../octTreeConcurrent/oFile
# cd ../octTreeConcurrent
# echo "3DPlummer"
# numactl -i all ./neighbors -d 3 -o oFile ../geometryData/data/3Dplummer_10M
# cd ../bench
# ./neighborsCheck -d 3 -k 1 -r 100 ../geometryData/data/3Dplummer_10M ../octTreeConcurrent/oFile
