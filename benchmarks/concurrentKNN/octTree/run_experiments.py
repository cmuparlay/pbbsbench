import argparse
import os
import sys
import multiprocessing

from create_graphs import *

# parameters:
#  - datastructures: neighbors_bench only for now
#  - threads, update_percent, number of nearest neighbors

ds_options = {
  "neighbors_bench" : "neighbors_bench",
  "neighbors_bench_path_copy" : "neighbors_bench_path_copy",
  "neighbors_bench_lockfree" : "neighbors_bench_lockfree",
  "neighbors_bench_hoh" : "neighbors_bench_hoh",
  "range_bench" : "../../rangeQueryKDTree/range/range_bench",
  "range_bench_path_copy" : "../../rangeQueryKDTree/range/range_bench_path_copy",
}

ds_keys = {
  "neighbors_bench" : "neighbors_bench",
  "neighbors_bench_path_copy" : "neighbors_bench_path_copy",
  "neighbors_bench_lockfree" : "neighbors_bench_lockfree",
  "neighbors_bench_hoh" : "neighbors_bench_hoh",
  "../../rangeQueryKDTree/range/range_bench" : "range_bench",
  "../../rangeQueryKDTree/range/range_bench_path_copy" : "range_bench_path_copy",
}


data_options = {
  "2DinCube_20M" : "../geometryData/data/2DinCube_20000000",
  "3DinCube_2M" : "../geometryData/data/3DinCube_2000000",
  "3DinCube_20M" : "../geometryData/data/3DinCube_20000000",
  "3Dplummer_20M" : "../geometryData/data/3Dplummer_20000000",
  "lucy3D_2M" : "/ssd0/angel/lucy3D_2M",
  "lucy3D_14M" : "/ssd0/angel/lucy3D_14M",
  "thai_statue2M" : "/ssd0/thai_statue/thai_statue2M",
  "thai_statue5M" : "/ssd0/thai_statue/thai_statue5M",
}

parser = argparse.ArgumentParser()
parser.add_argument("ds_type", help="datastructure type")
parser.add_argument("threads", help="Number of threads")
parser.add_argument("ratios", help="Update ratio, number between 0 and 100.")
parser.add_argument("query_sizes", help="size of query (k)")
parser.add_argument("dimension", help="dimension of data")
parser.add_argument("input_name", help="input name")
parser.add_argument("-t", "--test_only", help="test script",
                    action="store_true")
parser.add_argument("-g", "--graphs_only", help="graphs only",
                    action="store_true")
parser.add_argument("-p", "--paper_ver", help="paper version of graphs, no title or legends", action="store_true")

args = parser.parse_args()
print("datastructure: " + args.ds_type)
print("threads: " + args.threads)
print("update percent: " + args.ratios)
print("query_sizes: " + args.query_sizes)
print("dimension: " + args.dimension)
print("input_name: " + args.input_name)

test_only = args.test_only
graphs_only = args.graphs_only
rounds = 3

if test_only:
  rounds = 1

maxcpus = multiprocessing.cpu_count()
already_ran = set()

def string_to_list(s):
  s = s.strip().strip('[').strip(']').split(',')
  return [ss.strip() for ss in s]

def to_list(s):
  if type(s) == list:
    return s
  return [s]

def runstring(test, op, outfile, k):
    if op in already_ran:
        return
    already_ran.add(op)
    os.system("echo \"" + op + "\"")
    os.system("echo \"datastructure: " + ds_keys[test] + "-"+str(k)+"\"")
    os.system("echo \"" + op + "\" >> " + outfile)
    os.system("echo \"datastructure: " + ds_keys[test] + "-"+str(k) + "\" >> " + outfile)
    x = os.system(op + " >> " + outfile)
    if (x) :
        if (os.WEXITSTATUS(x) == 0) : raise NameError("  aborted: " + op)
        os.system("echo Failed")
    
def runtest(test,procs,u,k,d,infile,extra,outfile) :
    r = rounds
    otherargs = " -c -t 10.0 "

    runstring(test, "PARLAY_NUM_THREADS=" + str(min(int(procs), maxcpus)) + " numactl -i all ./" + test + " -r " + str(r) + " -d " + str(d) + " -k " + str(k) + " -p " + str(procs) + extra + " -u " + str(u) + otherargs + " " + infile, outfile, k)


exp_type = ""
ds_type = args.ds_type
threads = args.threads
ratios = args.ratios
query_sizes = args.query_sizes
dimension = args.dimension
input_name = args.input_name
input_file = data_options[input_name]

if '[' in args.threads:
  exp_type = "scalability"
  threads = string_to_list(threads)
else:
  print('invalid argument')
  exit(1)

if '[' in args.query_sizes:
  query_sizes = string_to_list(query_sizes)
else:
  print('invalid argument')
  exit(1)

if '[' in args.ds_type:
  ds_types = string_to_list(ds_type)
else:
  print('invalid argument')
  exit(1)


for i in ds_types:
  if i.find("range") != -1:
    exp_category = "range"
  else:
    exp_category = "neighbors"

outfile = "results/" + exp_category + "-".join([ratios+"up", input_name]) + ".txt"

datastructures=[]
for ds in ds_types:
  datastructures.append(ds_options[ds])

if not graphs_only:
  # clear output file
  os.system("echo \"\" > " + outfile)
  for ds in datastructures:
    for th in to_list(threads):
      for k in to_list(query_sizes):
        runtest(ds,th,ratios,k,dimension,input_file,"",outfile)

throughput = {}
stddev = {}
threads = []
ratios = []
algs = []

readResultsFile(outfile, throughput, stddev, threads, ratios, algs)

threads.sort()
ratios.sort()

print('threads: ' + str(threads))
print('update ratios: ' + str(ratios))
print('algs: ' + str(algs))
print(throughput)

graph_name=exp_category+input_name

alg_names=[]
for ds in ds_types:
  for k in query_sizes:
    alg_names.append(ds+"-"+str(k))

plot_scalability_graphs(throughput, stddev, threads, ratios, alg_names, graph_name, args.paper_ver)
