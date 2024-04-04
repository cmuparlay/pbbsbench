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
  "neighbors_bench_path_copy_lockfree" : "neighbors_bench_path_copy_lockfree",
  "neighbors_bench_lockfree" : "neighbors_bench_lockfree",
  "neighbors_bench_hoh" : "neighbors_bench_hoh",
  "working_set_bench" : "working_set_bench",
  "working_set_bench_lockfree" : "working_set_bench_lockfree",
  "working_set_bench_hoh" : "working_set_bench_hoh",
  "working_set_bench_path_copy" : "working_set_bench_path_copy",
  "range_bench" : "../../rangeQueryKDTree/range/range_bench",
  "range_bench_path_copy" : "../../rangeQueryKDTree/range/range_bench_path_copy",
  "LFKDTree" : "~/ConcurrentKDTree_o/dist/ConcurrentKDTree.jar",
}

ds_keys = {
  "neighbors_bench" : "neighbors_bench",
  "neighbors_bench_path_copy" : "neighbors_bench_path_copy",
  "neighbors_bench_path_copy_lockfree" : "neighbors_bench_path_copy_lockfree",
  "neighbors_bench_lockfree" : "neighbors_bench_lockfree",
  "neighbors_bench_hoh" : "neighbors_bench_hoh",
  "working_set_bench" : "working_set_bench",
  "working_set_bench_lockfree" : "working_set_bench_lockfree",
  "working_set_bench_hoh" : "working_set_bench_hoh",
  "working_set_bench_path_copy" : "working_set_bench_path_copy",
  "../../rangeQueryKDTree/range/range_bench" : "range_bench",
  "../../rangeQueryKDTree/range/range_bench_path_copy" : "range_bench_path_copy",
  "~/ConcurrentKDTree_o/dist/ConcurrentKDTree.jar" : "LFKDTree",
}


data_options = {
  "2DinCube_20M" : "../geometryData/data/2DinCube_20000000",
  "3DinCube_2K" : "../geometryData/data/3DinCube_2000",
  "3DinCube_20K" : "../geometryData/data/3DinCube_20000",
  "3DinCube_200K" : "../geometryData/data/3DinCube_200000",
  "3DinCube_2M" : "../geometryData/data/3DinCube_2000000",
  "3DinCube_20M" : "../geometryData/data/3DinCube_20000000",
  "3DinCube_WorkingSet_11M" : "../geometryData/data/3DinCube_WorkingSet_11M",
  "3DinCube_200M" : "/ssd0/geometryData/3DinCube_200M",
  "3DinCube_2B" : "/ssd0/geometryData/3DinCube_2B",
  "3Dplummer_20M" : "../geometryData/data/3Dplummer_20000000",
  "lucy3D_2M" : "/ssd0/angel/lucy3D_2M",
  "lucy3D_14M" : "/ssd0/angel/lucy3D_14M",
  "thai_statue2M" : "/ssd0/thai_statue/thai_statue2M",
  "thai_statue5M" : "/ssd0/thai_statue/thai_statue5M",
}

data_sizes = {
  "2DinCube_20M" : 10000000,
  "3DinCube_2K" : 1000,
  "3DinCube_20K" : 10000,
  "3DinCube_200K" : 100000,
  "3DinCube_2M" : 1000000,
  "3DinCube_20M" : 10000000,
  "3DinCube_200M" : 100000000,
  "3Dplummer_20M" : 10000000,
  "lucy3D_2M" : 1000000,
  "lucy3D_14M" : 7000000,
  "thai_statue2M" : 1000000,
  "thai_statue5M" : 2500000,
  "3DinCube_WorkingSet_11M" : 10000000,
}

filename_sizes = {
  "../geometryData/data/2DinCube_20000000" : 10000000,
  "../geometryData/data/3DinCube_2000" : 1000,
  "../geometryData/data/3DinCube_20000" : 10000,
  "../geometryData/data/3DinCube_200000" : 100000,
  "../geometryData/data/3DinCube_2000000": 1000000,
  "../geometryData/data/3DinCube_20000000" : 10000000,
  "../geometryData/data/3DinCube_WorkingSet_11M" : 10000000,
  "/ssd0/geometryData/3DinCube_200M" : 100000000,
  "/ssd0/geometryData/3DinCube_2B" : 1000000000,
  "../geometryData/data/3Dplummer_20000000" : 10000000,
  "/ssd0/angel/lucy3D_2M" : 1000000,
  "/ssd0/angel/lucy3D_14M" : 7000000,
  "/ssd0/thai_statue/thai_statue2M" : 1000000,
  "/ssd0/thai_statue/thai_statue5M" : 2500000,
}




parser = argparse.ArgumentParser()
parser.add_argument("ds_type", help="datastructure type")
parser.add_argument("threads", help="Number of threads")
parser.add_argument("ratios", help="Update ratio, number between 0 and 100.")
parser.add_argument("query_sizes", help="size of query (k)")
parser.add_argument("dimension", help="dimension of data")
parser.add_argument("input_names", help="input names")
parser.add_argument("-t", "--test_only", help="test script",
                    action="store_true")
parser.add_argument("-g", "--graphs_only", help="graphs only",
                    action="store_true")
parser.add_argument("-p", "--paper_ver", help="paper version of graphs, no title or legends", action="store_true")

args = parser.parse_args()
print("datastructure: " + args.ds_type)
print("threads: " + args.threads)
print("update percents: " + args.ratios)
print("query_sizes: " + args.query_sizes)
print("dimension: " + args.dimension)
print("input_names: " + args.input_names)

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
    if(int(procs) < maxcpus): 
      tr = maxcpus
    else:
      tr = int(procs)

    runstring(test, "PARLAY_NUM_THREADS=" + str(tr) + " numactl -i all ./" + test + " -r " + str(r) + " -d " + str(d) + " -k " + str(k) + " -p " + str(procs) + extra + " -u " + str(u) + otherargs + " " + infile, outfile, k)

def runCompetitor(test,procs,u,d,n,outfile) :
    s = 100-int(u)
    op = "java -jar " + test + " -a LFKDTree -d " + str(d) + " -n " + str(procs) + " -r " + str(s) + " -i " + str(int(int(u)/2)) + " -x " + str(int(int(u)/2)) + " -k "+ str(n) + " -m " + str(d)
    os.system("echo \"" + op + "\"")
    os.system("echo \"datastructure: " + ds_keys[test] + "-"+str(k)+"\"")
    os.system("echo \"" + op + "\" >> " + outfile)
    os.system("echo \"datastructure: " + ds_keys[test] + "-"+str(k) + "\" >> " + outfile)
    os.system("echo \"" + op + "\" >> " + outfile)
    x = os.system(op + " >> " + outfile)
    if (x) :
        if (os.WEXITSTATUS(x) == 0) : raise NameError("  aborted: " + op)
        os.system("echo Failed")
    


exp_type = ""
ds_type = args.ds_type
threads = args.threads
ratios = args.ratios
query_sizes = args.query_sizes
dimension = args.dimension
input_names = args.input_names

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

if '[' in args.ratios:
  ratios = string_to_list(ratios)
else:
  print('invalid argument')
  exit(1)

if '[' in args.input_names:
  input_names = string_to_list(input_names)
else:
  print('invalid argument')
  exit(1)

input_files = []
for name in input_names:
  input_files.append(data_options[name])
print(input_files)



for i in ds_types:
  if i.find("range") != -1:
    exp_category = "range"
  else:
    exp_category = "neighbors"

outfile = "results/" + exp_category + "-".join([input_names[0]]) + ".txt"

datastructures=[]
for ds in ds_types:
  datastructures.append(ds_options[ds])

if not graphs_only:
  # clear output file
  os.system("echo \"\" > " + outfile)
  for file in input_files:
    for ds in datastructures:
      for th in to_list(threads):
        for k in to_list(query_sizes):
          for u in to_list(ratios):
            if((ds_keys[ds] == "LFKDTree") and (int(k) == 1)):
              runCompetitor(ds,th,u,dimension,filename_sizes[file],outfile)
            elif((ds_keys[ds] == "LFKDTree")):
              continue
            else:
              runtest(ds,th,u,k,dimension,file,"",outfile)

throughput = {}
stddev = {}
threads = []
ratios = []
algs = []
sizes = []

readResultsFile(outfile, throughput, stddev, threads, ratios, sizes, algs)

threads.sort()
ratios.sort()

print('threads: ' + str(threads))
print('update ratios: ' + str(ratios))
print('algs: ' + str(algs))

graph_name=exp_category+input_names[0]

alg_names=[]
for ds in ds_types:
  for k in query_sizes:
    alg_names.append(ds+"-"+str(k))

args.paper_ver = True
if(len(input_names) <= 1):
  plot_scalability_graphs(throughput, stddev, threads, ratios, sizes[0], alg_names, "scalability_"+graph_name, args.paper_ver)
  plot_ratio_graphs(throughput, stddev, threads, ratios, sizes[0], alg_names, "ratio_"+graph_name, args.paper_ver)
else:
  plot_size_graphs(throughput, stddev, threads, ratios, sizes, alg_names, "sizes_"+graph_name, args.paper_ver)


