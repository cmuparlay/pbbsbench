#!/bin/bash
cd ~/big-ann-benchmarks
# # download datasets
# #may need to be run a couple times, occasionally the connection returns an error
python3 create_dataset.py --dataset bigann-10M 
python3 create_dataset.py --dataset bigann-100M 
python3 create_dataset.py --dataset bigann-1B 
python3 create_dataset.py --dataset msspacev-1M 
python3 create_dataset.py --dataset msspacev-10M 
python3 create_dataset.py --dataset msspacev-100M 
python3 create_dataset.py --dataset msspacev-1B
python3 create_dataset.py --dataset text2image-1M 
python3 create_dataset.py --dataset text2image-10M 
python3 create_dataset.py --dataset text2image-100M
python3 create_dataset.py --dataset text2image-1B 
python3 create_dataset.py --dataset ssnpp-1M 
python3 create_dataset.py --dataset ssnpp-10M 
python3 create_dataset.py --dataset ssnpp-100M
python3 create_dataset.py --dataset ssnpp-1B
