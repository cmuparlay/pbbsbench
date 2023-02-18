cd ~/big-ann-benchmarks
P=~/pbbsbench/benchmarks/ANN/scripts/FAISS
datasets=("bigann-10M")
for d in ${datasets[@]}; do
  python3 run.py --definitions $P/ANN.yaml --algorithm faiss-t1 --dataset $d
done
sudo chmod -R 777 results/
python data_export.py --output res.csv