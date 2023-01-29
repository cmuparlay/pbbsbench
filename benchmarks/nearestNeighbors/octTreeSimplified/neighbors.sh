#!/bin/bash

# in checkNeighbors, r = number of points to calculate exhaustive nbh for

make 

echo "2DinCube"
numactl -i all ./neighbors -d 2 -k 1 -o oFile ../geometryData/data/2DinCube_10M
cd ../bench
./neighborsCheck -d 2 -k 1 -r 100 ../geometryData/data/2DinCube_10M ../octTreeSimplified/oFile
cd ../octTreeSimplified
echo "2DKuzmin"
numactl -i all ./neighbors -d 2 ../geometryData/data/2Dkuzmin_10M
cd ../bench
./neighborsCheck -d 2 -k 1 -r 100 ../geometryData/data/2Dkuzmin_10M ../octTreeSimplified/oFile
cd ../octTreeSimplified
echo "3DonSphere"
numactl -i all ./neighbors -d 3 ../geometryData/data/3DonSphere_10M
cd ../bench
./neighborsCheck -d 3 -k 1 -r 100 ../geometryData/data/3DonSphere_10M ../octTreeSimplified/oFile
cd ../octTreeSimplified
echo "3DinCube"
numactl -i all ./neighbors -d 3 ../geometryData/data/3DinCube_10M
cd ../bench
./neighborsCheck -d 3 -k 1 -r 100 ../geometryData/data/3DinCube_10M ../octTreeSimplified/oFile
cd ../octTreeSimplified
echo "3DPlummer"
numactl -i all ./neighbors -d 3 ../geometryData/data/3Dplummer_10M
cd ../bench
./neighborsCheck -d 3 -k 1 -r 100 ../geometryData/data/3Dplummer_10M ../octTreeSimplified/oFile
cd ../octTreeSimplified