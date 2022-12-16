#!/bin/bash

make clean all

echo "2DinCube"
echo ''
echo ''
./neighbors -d 2 ../geometryData/data/2DinCube_10M
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