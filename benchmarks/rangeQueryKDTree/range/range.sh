#!/bin/bash

# in checkNeighbors, r = number of points to calculate exhaustive nbh for

make 

echo "2DinCube"
./range -d 2 -rad .0008 -o oFile ../geometryData/data/2DinCube_10M
cd ../bench
./rangeCheck -d 2 -rad .0008 -r 100 ../geometryData/data/2DinCube_10M ../range/oFile
cd ../range
# echo "2DKuzmin"
numactl -i all ./range -d 2 -rad .002 -o oFile ../geometryData/data/2Dkuzmin_10M
cd ../bench
./rangeCheck -d 2 -rad .002 -r 100 ../geometryData/data/2Dkuzmin_10M ../range/oFile
cd ../range
echo "3DonSphere"
numactl -i all ./range -d 3 -rad .0014 -o oFile ../geometryData/data/3DonSphere_10M
cd ../bench
./rangeCheck -d 3 -rad .0014 -r 100 ../geometryData/data/3DonSphere_10M ../range/oFile
cd ../range
echo "3DinCube"
numactl -i all ./range -d 3 -rad .01 -o oFile ../geometryData/data/3DinCube_10M
cd ../bench
./rangeCheck -d 3 -rad .01 -r 100 ../geometryData/data/3DinCube_10M ../range/oFile
cd ../range
echo "3DPlummer"
numactl -i all ./range -d 3 -rad .014 -o oFile ../geometryData/data/3Dplummer_10M
cd ../bench
./rangeCheck -d 3 -rad .014 -r 100 ../geometryData/data/3Dplummer_10M ../range/oFile