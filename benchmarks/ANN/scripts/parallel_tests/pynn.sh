#/bin/bash

for i in 1 2 8 12 24 48 96 144 192; do 
    export NUMBA_NUM_THREADS=$i
    python3 pynn.py >> pynn_data.txt
done