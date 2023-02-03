Here I describe how to use the scripts in this repository.

PREREQUISITES:
This assumes that you have the repos for pbbsbench and big-ann-benchmarks pulled from github.
For big-ann-benchmarks, use [these instructions](https://github.com/harsha-simhadri/big-ann-benchmarks).
For pbbsbench, use [these instructions](https://cmuparlay.github.io/pbbsbench/).
It assumes that you have an SSD mounted and named /ssd1 (at least 4 TB needed). In big-ann-benchmarks, the "data" folder should be a symlink to /ssd1/data. You should also create a directory /ssd1/results and create a directory with the name of each dataset (bigann, MSSPACEV1B, text2image1B, FB_ssnpp). This is where the CSV file results and graphs wil be stored. 

Now that you've done all this, you can get started running scripts in this repo. 

First, run "download_datasets.sh." It is isolated in its own script because in my experience the downloads sometimes fail. It may need to be re-ran two or three times.  This needs 2TB for data prep, 4TB for all indices.  

Now there are four more scripts contained in "all_experiments.sh." 

First, "prepare_datasets.sh" computes groundtruth files for all the 1M slices, since they are not contained in big-ann-benchmarks. It also computes the bigann-1M slice. 

Three files of experiments follow: "vamana.sh," "hcnng.sh," and "pynndescent.sh." Each experiment builds a number of ANN graphs (6 billion scale builds per script, with the exception of pynndescent which only goes to 100 million) graph, queries it, and then writes the graph to disk in /ssd1/results. It also writes a CSV of statistics (QPS, recall, etc.) to /ssd1/results. If you don't want to store the graph, you can remove this option by removing the path in "-o" in each command. 

