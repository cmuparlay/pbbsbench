#!/bin/bash

make 

echo "2DinCube"
echo ''
echo ''
./neighbors -d 2 -k 1 -o oFile ../geometryData/data/2DinCube_10M
cd ../bench
./neighborsCheck -d 2 -k 1 -r 100 ../geometryData/data/2DinCube_10M ../octTreeSimplified/oFile
# echo "2Dkuzmin"
# echo ''
# echo ''
# ./neighbors -d 2 ../geometryData/data/2Dkuzmin_10M
# echo "3DonSphere"
# echo ''
# echo ''
# ./neighbors -d 3 ../geometryData/data/3DonSphere_10M
# echo "3DinCube"
# echo ''
# echo ''
# ./neighbors -d 3 ../geometryData/data/3DinCube_10M
# echo "3DPlummer"
# echo ''
# echo ''
# ./neighbors -d 3 ../geometryData/data/3Dplummer_10M