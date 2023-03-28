#!/bin/bash
cd ~/big-ann-benchmarks
rm results
mkdir /ssd1/results/FAISSresults3
ln -s /ssd1/results/FAISSresults3 results
P=~/pbbsbench/benchmarks/ANN/scripts/FAISS
R=/ssd1/results/FAISSresults3

# nohup python3 run.py --definitions $P/ANN2.yaml --algorithm faiss-t1 --dataset "bigann-10M" > $R/bigann-10M2.log
# nohup python3 run.py --definitions $P/ANN2.yaml --algorithm faiss-t1 --dataset "msspacev-10M" > $R/msspacev-10M2.log
# nohup python3 run.py --definitions $P/ANN2.yaml --algorithm faiss-t1 --dataset "text2image-10M" > $R/text2image-10M2.log
# nohup python3 run.py --definitions $P/ANN2.yaml --algorithm faiss-t1 --dataset "ssnpp-10M" > $R/ssnpp-10M2.log

# nohup python3 run.py --definitions $P/ANN2.yaml --algorithm faiss-t1 --dataset "bigann-100M" > $R/bigann-100M2.log
nohup python3 run.py --definitions $P/ANN2.yaml --algorithm faiss-t1 --dataset "msspacev-100M" --rebuild > $R/msspacev-100M2.log
# nohup python3 run.py --definitions $P/ANN2.yaml --algorithm faiss-t1 --dataset "text2image-100M" > $R/text2image-100M2.log
# nohup python3 run.py --definitions $P/ANN2.yaml --algorithm faiss-t1 --dataset "ssnpp-100M" > $R/ssnpp-100M2.log

# nohup python3 run.py --definitions $P/ANN2.yaml --algorithm faiss-t1 --dataset "bigann-1B" > $R/bigann-1B2.log
nohup python3 run.py --definitions $P/ANN2.yaml --algorithm faiss-t1 --dataset "msspacev-1B" --rebuild > $R/msspacev-1B2.log
# nohup python3 run.py --definitions $P/ANN2.yaml --algorithm faiss-t1 --dataset "text2image-1B" > $R/text2image-1B2.log
# nohup python3 run.py --definitions $P/ANN2.yaml --algorithm faiss-t1 --dataset "ssnpp-1B" > $R/ssnpp-1B2.log

sudo chmod -R 777 results/
python data_export.py --output /ssd1/results/FAISSresults3/res3.csv