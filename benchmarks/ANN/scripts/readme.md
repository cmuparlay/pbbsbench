Here I describe how to use the scripts in this repository.

PREREQUISITES:
This assumes that you have the repos for pbbsbench and big-ann-benchmarks pulled from github.
For big-ann-benchmarks, use [these instructions](https://github.com/harsha-simhadri/big-ann-benchmarks).
For pbbsbench, use [these instructions](https://cmuparlay.github.io/pbbsbench/).
It assumes that you have an SSD mounted and named /ssd1 (at least 4 TB needed). In big-ann-benchmarks, the "data" folder should be a symlink to /ssd1/data. You should also create a directory /ssd1/results where CSV results and graphs will be stored, or you can modify the paths in the script as needed.

Now that you've done all this, you can get started running scripts in this repo. 

First, run "download_datasets.sh." It is isolated in its own script because in my experience the downloads sometimes fail. It may need to be re-ran two or three times. 

Now there are four more scripts contained in "all_experiments.sh." 

First, "prepare_datasets.sh" computes groundtruth files for all the 1M slices, since they are not contained in big-ann-benchmarks. It also computes the bigann-1M slice. 

Three files of experiments follow: "vamana.sh," "hcnng.sh," and "pynndescent.sh." Each experiment builds a graph, queries it, and then writes the graph to disk in /ssd1/results. It also writes a CSV of statistics (QPS, recall, etc.) to /ssd1/results. If you don't want to store the graph, you can remove this option by removing the path in "-o" in each command. 

