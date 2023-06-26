import pynndescent
import numpy as np
import h5py
from urllib.request import urlretrieve
import os
import time

def get_ann_benchmark_data(dataset_name):
    if not os.path.exists(f"{dataset_name}.hdf5"):
        print(f"Dataset {dataset_name} is not cached; downloading now ...")
        urlretrieve(f"http://ann-benchmarks.com/{dataset_name}.hdf5", f"{dataset_name}.hdf5")
    hdf5_file = h5py.File(f"{dataset_name}.hdf5", "r")
    return np.array(hdf5_file['train']), np.array(hdf5_file['test']), hdf5_file.attrs['distance']

sift_train, sift_test, _ = get_ann_benchmark_data('sift-128-euclidean')
sift_train.shape
sift_test.shape

print("Got data, beginning index build")
tic = time.perf_counter()
index = pynndescent.NNDescent(sift_train, n_neighbors=60, diversify_prob=.8, pruning_degree_multiplier=2.0)
# index = pynndescent.NNDescent(sift_train)
toc = time.perf_counter()
print(f"Built index in {toc - tic:0.4f} seconds")