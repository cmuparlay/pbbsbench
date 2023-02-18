#!/bin/bash
cd ~/big-ann-benchmarks
P=~/pbbsbench/benchmarks/ANN/scripts/FAISS

python3 run.py --definitions $P/ANN.yaml --algorithm faiss-t1 --dataset "bigann-10M"
python3 run.py --definitions $P/ANN.yaml --algorithm faiss-t1 --dataset "msspacev-10M"
python3 run.py --definitions $P/ANN.yaml --algorithm faiss-t1 --dataset "text2image-10M"
python3 run.py --definitions $P/ANN.yaml --algorithm faiss-t1 --dataset "ssnpp-10M"

python3 run.py --definitions $P/ANN.yaml --algorithm faiss-t1 --dataset "bigann-100M"
python3 run.py --definitions $P/ANN.yaml --algorithm faiss-t1 --dataset "msspacev-100M"
python3 run.py --definitions $P/ANN.yaml --algorithm faiss-t1 --dataset "text2image-100M"
python3 run.py --definitions $P/ANN.yaml --algorithm faiss-t1 --dataset "ssnpp-100M"

python3 run.py --definitions $P/ANN.yaml --algorithm faiss-t1 --dataset "bigann-1B"
python3 run.py --definitions $P/ANN.yaml --algorithm faiss-t1 --dataset "msspacev-1B"
python3 run.py --definitions $P/ANN.yaml --algorithm faiss-t1 --dataset "text2image-1B"
python3 run.py --definitions $P/ANN.yaml --algorithm faiss-t1 --dataset "ssnpp-1B"

sudo chmod -R 777 results/
python data_export.py --output res.csv