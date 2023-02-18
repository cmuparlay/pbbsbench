cd ~/big-ann-benchmarks
P=~/pbbsbench/benchmarks/ANN/scripts/FAISS
datasets=("bigann-10M" "ssnpp-10M" "msspacev-10M" "text2image-10M" "bigann-100M" "ssnpp-100M" "msspacev-100M" "text2image-100M" "bigann-1B" "ssnpp-1B" "msspacev-1B" "text2image-1B")
for d in ${datasets[@]}; do
  python3 run.py --definitions $P/ANN.yaml --algorithm faiss-t1 --dataset $d
done
sudo chmod -R 777 results/
python data_export.py --output res.csv