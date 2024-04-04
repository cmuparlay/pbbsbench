import csv
import matplotlib as mpl
# mpl.use('Agg')
mpl.rcParams['grid.linestyle'] = ":"
mpl.rcParams.update({'font.size': 20})
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import os
import statistics as st

# paper_ver = False

class DSInfo:
    def __init__(self, color, marker, linestyle, name, binary, ds_type):
        self.color = color
        self.marker = marker
        self.linestyle = linestyle
        self.name = name
        self.ds_type = ds_type
        self.binary = binary

mk = ['o', 'v', '^', '1', 's', '+', 'x', 'D', '|', '>', '<',]

dsinfo = {'neighbors_bench-1' :  DSInfo("C0", mk[0], "-", "CLEANN-Tree, k=1", "blah", "blah"),
          'neighbors_bench-10' : DSInfo("C1", mk[1], "-", "CLEANN-Tree, k=10", "blah", "blah"),
          'neighbors_bench_hoh-1' :  DSInfo("C5", mk[5], "-", "CLEANN-Tree-HOH, k=1", "blah", "blah"),
          'neighbors_bench_hoh-10' : DSInfo("C6", mk[6], "-", "CLEANN-Tree-HOH, k=10", "blah", "blah"),
          'working_set_bench-1' :  DSInfo("C0", mk[0], "-", "CLEANN-Tree, k=1", "blah", "blah"),
          'working_set_bench-10' : DSInfo("C1", mk[1], "-", "CLEANN-Tree, k=10", "blah", "blah"),
          'working_set_bench_lockfree-1' :  DSInfo("C7", mk[7], "-", "CLEANN-Tree-LF, k=1", "blah", "blah"),
          'working_set_bench_lockfree-10' : DSInfo("C8", mk[8], "-", "CLEANN-Tree-LF, k=10", "blah", "blah"),
          'working_set_bench_hoh-1' :  DSInfo("C5", mk[5], "-", "CLEANN-Tree-HOH, k=1", "blah", "blah"),
          'working_set_bench_hoh-10' : DSInfo("C6", mk[6], "-", "CLEANN-Tree-HOH, k=10", "blah", "blah"),
          'working_set_bench_path_copy-1' :  DSInfo("C3", mk[3], "-", "CLEANN-Tree-PC, k=1", "blah", "blah"),
          'working_set_bench_path_copy-10' : DSInfo("C4", mk[4], "-", "CLEANN-Tree-PC, k=10", "blah", "blah"),
          'neighbors_bench_path_copy-1' :  DSInfo("C3", mk[3], "-", "CLEANN-Tree-PC, k=1", "blah", "blah"),
          'neighbors_bench_path_copy-10' : DSInfo("C4", mk[4], "-", "CLEANN-Tree-PC, k=10", "blah", "blah"),
          'neighbors_bench_lockfree-1' :  DSInfo("C7", mk[7], "-", "CLEANN-Tree-LF, k=1", "blah", "blah"),
          'neighbors_bench_lockfree-10' : DSInfo("C8", mk[8], "-", "CLEANN-Tree-LF, k=10", "blah", "blah"),
          'neighbors_bench_path_copy_lockfree-1' :  DSInfo("C2", mk[2], "-", "CLEANN-Tree-PC-LF, k=1", "blah", "blah"),
          'neighbors_bench_path_copy_lockfree-10' : DSInfo("C9", mk[9], "-", "CLEANN-Tree-PC-LF, k=10", "blah", "blah"),
          'range_bench-.014' : DSInfo("C0", mk[0], "-", "CLEANN-Tree, r=5", "blah", "blah"),
          'range_bench-.0176' : DSInfo("C1", mk[1], "-", "CLEANN-Tree, r=10", "blah", "blah"),
          'range_bench_path_copy-.014' : DSInfo("C3", mk[3], "-", "CLEANN-Tree-PC, r=5", "blah", "blah"),
          'range_bench_path_copy-.0176' : DSInfo("C4", mk[4], "-", "CLEANN-Tree-PC, r=10", "blah", "blah"),
          'LFKDTree-1' :  DSInfo("C2", mk[9], "-", "LFKDTree, k=1", "blah", "blah"),
} 

ratiomarkers = {0 : mk[0], 50 : mk[1], 100 : mk[2]}


def toString(algname, th, ratio, size):
    return algname + '-' + str(th) + 'u-' + str(ratio) + "s-" + str(size)


def export_legend(legend, filename="legend.pdf"):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(filename, dpi="figure", bbox_inches=bbox)

def avg(numlist):
  # if len(numlist) == 1:
  #   return numlist[0]
  total = 0.0
  length = 0
  for num in numlist:
    length=length+1
    total += float(num)
  if length > 0:
    return 1.0*total/length
  else:
    return -1

def readResultsFile(filename, throughput, stddev, threads, ratios, sizes, algs):
  throughputRaw = {}
  alg = ""
  th = ""
  ratio = ""
  maxkey = ""
  alpha = ""
  size = ""
  warm_up_runs = 0

  file = open(filename, 'r')
  query_sizes=[]
  for line in file.readlines():
    if line.find('query size:') != -1:
      k = int(line.split(': ')[1])
      query_sizes.append(k)

  # read csv into throughputRaw
  file = open(filename, 'r')
  for line in file.readlines():
    line = line.strip()
    if line.find('warm_up_runs') != -1:
      warm_up_runs = int(line.split(' ')[1])
    elif line.find('datastructure:') != -1:
      alg = line.split(' ')[1]
    elif line.find('threads:') != -1:
      th = int(line.split(' ')[1])
    elif line.find('update_percent:') != -1:
      ratio = int(line.split(' ')[1])
    elif line.find('index_size:') != -1:
      size = int(line.split(' ')[1])
    elif line.find('throughput (Mop/s):') != -1:
      tp = float(line.split(': ')[1])
      if alg not in algs:
        algs.append(alg)
      if th not in threads:
        threads.append(th)
      if ratio not in ratios:
        ratios.append(ratio)
      if size not in sizes:
        sizes.append(size)
      key = toString(alg, th, ratio, size)
      if key not in throughputRaw:
        throughputRaw[key] = []
      throughputRaw[key].append(tp)

  print(throughputRaw)
  # Average througputRaw into throughput

  for key in throughputRaw:
    results = throughputRaw[key][warm_up_runs:]
    throughput[key] = avg(results)
    stddev[key] = st.pstdev(results)

#   print(throughput)

# def plot_alpha_graph(throughput, stddev, thread, ratio, maxkey, alphas, algs, graph_name, paper_ver=False):
#   # print(graphtitle)
#   graphtitle = graph_name + '-' + str(thread) + 'th-' + str(maxkey) + 'size-' + str(ratio) + 'up'
#   if paper_ver:
#     outputFile = 'graphs/' + graphtitle.replace('.', '') + '.pdf'
#     mpl.rcParams.update({'font.size': 25})
#   else:
#     outputFile = 'graphs/' + graphtitle + '.png'

#   ymax = 0
#   series = {}
#   error = {}
#   for alg in algs:
#     # if (alg == 'BatchBST64' or alg == 'ChromaticBatchBST64') and (bench.find('-0rq') == -1 or bench.find('2000000000') != -1):
#     #   continue
#     # if toString3(alg, 1, bench) not in results:
#     #   continue
#     series[alg] = []
#     error[alg] = []
#     for alpha in alphas:
#       key = toString(alg, thread, ratio, maxkey, alpha)
#       if key not in throughput:
#         del series[alg]
#         del error[alg]
#         break
#       series[alg].append(throughput[key])
#       error[alg].append(stddev[key])
  
#   if len(series) < 3:
#     return

#   fig, axs = plt.subplots()
#   # fig = plt.figure()
#   opacity = 0.8
#   rects = {}
  
#   xpos = np.arange(start=0, stop=len(alphas))

#   for alg in algs:
#     alginfo = dsinfo[alg]
#     if alg not in series:
#       continue
#     ymax = max(ymax, max(series[alg]))
#     # rects[alg] = axs.errorbar(xpos, series[alg], yerr=error[alg],
#     # if graph_name == 'try_vs_strict_lock' and alginfo.name.find('leaftree-strict') == -1:
#     if 'leaftree-lf' in algs and alginfo.name.find('leaftree-strict') == -1:
#       rects[alg] = axs.plot(xpos, series[alg],
#       alpha=opacity,
#       color=alginfo.color,
#       linewidth=3.0,
#       #hatch=hatch[ds],
#       linestyle=alginfo.linestyle,
#       marker=alginfo.marker,
#       markersize=14,
#       label=alginfo.name.replace('leaftree-', 'leaftree-trylock-'))
#     else:
#       rects[alg] = axs.plot(xpos, series[alg],
#         alpha=opacity,
#         color=alginfo.color,
#         linewidth=3.0,
#         #hatch=hatch[ds],
#         linestyle=alginfo.linestyle,
#         marker=alginfo.marker,
#         markersize=14,
#         label=alginfo.name)

#   plt.xticks(xpos, alphas)
#   axs.set_ylim(bottom=-0.02*ymax)
#   # plt.xticks(threads, threads)
#   # axs.set_xlabel('Number of threads')
#   # axs.set_ylabel('Throughput (Mop/s)')
#   axs.set(xlabel='Zipfian parameter', ylabel='Throughput (Mop/s)')
#   legend_x = 1
#   legend_y = 0.5 
#   # if this_file_name == 'Update_heavy_with_RQ_-_100K_Keys':
#   #   plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))

#   # plt.legend(framealpha=0.0)

#   plt.grid()
#   if 'leaftree-lf' in algs:
#     plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))
#   elif not paper_ver:
#     plt.title(graphtitle)
#     plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))
#   plt.savefig(outputFile, bbox_inches='tight')
#   plt.close('all')

# def plot_alpha_graphs(throughput, stddev, threads, ratios, maxkeys, alphas, algs, graph_name, paper_ver=False):
#   for thread in threads:
#     for size in maxkeys:
#       for ratio in ratios:
#         sufficient_datapoints = False
#         for alg in algs:
#           num_datapoints = 0
#           for alpha in alphas:
#             if toString(alg, thread, ratio, size, alpha) in throughput:
#               num_datapoints += 1
#           if num_datapoints > 2:
#             sufficient_datapoints = True
#         if sufficient_datapoints:
#           plot_alpha_graph(throughput, stddev, thread, ratio, size, alphas, algs, graph_name, paper_ver)

def plot_size_graph(throughput, stddev, thread, ratios, maxkeys, algs, graph_name, paper_ver=False):
  # print(graphtitle)
  graphtitle = graph_name + '-' + str(thread) + 'th'
  if paper_ver:
    outputFile = 'graphs/' + graphtitle.replace('.', '') + '.pdf'
    mpl.rcParams.update({'font.size': 25})
  else:
    outputFile = 'graphs/' + graphtitle + '.png'

  ymax = 0
  series = {}
  error = {}
  for alg in algs:
    for ratio in ratios:
      series[alg] = []
      error[alg] = []
      for maxkey in maxkeys:
        print(maxkey)
        key = toString(alg, thread, ratio, maxkey)
        if key not in throughput:
          print("breaking")
          del series[alg]
          del error[alg]
          break
        series[alg].append(throughput[key])
        error[alg].append(stddev[key])
  # if len(series) < 3:
  #   print("series too short; returning")
  #   return

  fig, axs = plt.subplots()
  # fig = plt.figure()
  opacity = 0.8
  rects = {}
  
  for alg in algs:
    for ratio in ratios:
      alginfo = dsinfo[alg]
      if alg not in series:
        continue
      ymax = max(ymax, max(series[alg]))
      # rects[alg] = axs.errorbar(maxkeys, series[alg], yerr=error[alg],
      rects[alg] = axs.plot(maxkeys, series[alg],
        alpha=opacity,
        color=alginfo.color,
        #hatch=hatch[ds],
        linewidth=3.0,
        linestyle=alginfo.linestyle,
        marker=ratiomarkers[ratio],
        markersize=14,
        label=alginfo.name+",u="+str(ratio))

  # if maxkeys[-1] > 1000000:
  #   plt.axvline(1000000, linestyle='--', color='grey') 
  axs.set_xscale('log')
  axs.set_ylim(bottom=-0.02*ymax)
  # plt.xticks(threads, threads)
  # axs.set_xlabel('Number of threads')
  # axs.set_ylabel('Throughput (Mop/s)')
  axs.set(xlabel='Datastructure size', ylabel='Throughput (Mop/s)')
  legend_x = 1
  legend_y = 0.5 
  # if this_file_name == 'Update_heavy_with_RQ_-_100K_Keys':
  #   plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))

  # plt.legend(framealpha=0.0)
  plt.grid()
  if not paper_ver:
    plt.title(graphtitle)
    plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))
  plt.savefig(outputFile, bbox_inches='tight')

  if paper_ver:
    legend = plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y), ncol=8, framealpha=0.0)
    export_legend(legend, 'graphs/' + graph_name + '_legend.pdf')
  plt.close('all')

def plot_size_graphs(throughput, stddev, threads, ratios, maxkeys, algs, graph_name, paper_ver=False):
  for thread in threads:
    sufficient_datapoints = False
    xaxis = []
    for alg in algs:
      for ratio in ratios:
        num_datapoints = 0
        for size in maxkeys:
          if toString(alg, thread, ratio, size) in throughput:
            num_datapoints += 1
            if size not in xaxis:
              xaxis.append(size)
        if num_datapoints > 2:
          sufficient_datapoints = True
    if sufficient_datapoints:
      print("Sufficient datapoints")
      plot_size_graph(throughput, stddev, thread, ratios, xaxis, algs, graph_name, paper_ver)


def plot_ratio_graph(throughput, stddev, thread, ratios, size, algs, graph_name, paper_ver=False):
  # print(graphtitle)
  graphtitle = graph_name + '-' + str(thread) + 'th'
  if paper_ver:
    outputFile = 'graphs/' + graphtitle.replace('.', '') + '.pdf'
    mpl.rcParams.update({'font.size': 25})
  else:
    outputFile = 'graphs/' + graphtitle + '.png'
  print(outputFile)

  ymax = 0
  series = {}
  error = {}
  for alg in algs:
    series[alg] = []
    error[alg] = []
    for ratio in ratios:
      key = toString(alg, thread, ratio, size)
      if key not in throughput:
        del series[alg]
        del error[alg]
        print("breaking")
        break
      series[alg].append(throughput[key])
      error[alg].append(stddev[key])

  fig, axs = plt.subplots()
  # fig = plt.figure()
  opacity = 0.8
  rects = {}
  
  xpos = np.arange(start=0, stop=len(ratios))

  for alg in algs:
    alginfo = dsinfo[alg]
    if alg not in series:
      continue
    ymax = max(ymax, max(series[alg]))
    # rects[alg] = axs.errorbar(xpos, series[alg], yerr=error[alg],
    rects[alg] = axs.plot(xpos, series[alg],
      alpha=opacity,
      color=alginfo.color,
      linewidth=3.0,
      #hatch=hatch[ds],
      linestyle=alginfo.linestyle,
      marker=alginfo.marker,
      markersize=14,
      label=alginfo.name)

  plt.xticks(xpos, ratios)
  axs.set_ylim(bottom=-0.02*ymax)
  # plt.xticks(threads, threads)
  # axs.set_xlabel('Number of threads')
  # axs.set_ylabel('Throughput (Mop/s)')
  axs.set(xlabel='Update percentage', ylabel='Throughput (Mop/s)')
  legend_x = 1
  legend_y = 0.5 
  # if this_file_name == 'Update_heavy_with_RQ_-_100K_Keys':
  #   plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))

  # plt.legend(framealpha=0.0)
  plt.grid()
  if not paper_ver:
    plt.title(graphtitle)
    plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))
  
  plt.savefig(outputFile, bbox_inches='tight')

  if paper_ver:
    legend = plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y), ncol=8, framealpha=0.0)
    export_legend(legend, 'graphs/' + graph_name + '_legend.pdf')
  plt.close('all')

def plot_ratio_graphs(throughput, stddev, threads, ratios, size, algs, graph_name, paper_ver=False):
  for thread in threads:
    sufficient_datapoints = False
    for alg in algs:
      num_datapoints = 0
      for ratio in ratios:
        if toString(alg, thread, ratio, size) in throughput:
          num_datapoints += 1
      if num_datapoints > 2:
        sufficient_datapoints = True
    if sufficient_datapoints:
      plot_ratio_graph(throughput, stddev, thread, ratios, size, algs, graph_name, paper_ver)


def plot_scalability_graph(throughput, stddev, threads, ratio, size, algs, graph_name, paper_ver=False):
  graphtitle = graph_name + '-' + str(ratio) 
  if paper_ver:
    outputFile = 'graphs/' + graphtitle.replace('.', '') + '.pdf'
    mpl.rcParams.update({'font.size': 25})
  else:
    outputFile = 'graphs/' + graphtitle + '.png'
  print("plotting " + outputFile)

  ymax = 0
  series = {}
  error = {}
  for alg in algs:
    series[alg] = []
    error[alg] = []
    for th in threads:
      key = toString(alg, th, ratio, size)
      if key not in throughput:
        del series[alg]
        del error[alg]
        break
      series[alg].append(throughput[key])
      error[alg].append(stddev[key])
  fig, axs = plt.subplots()
  # fig = plt.figure()
  opacity = 0.8
  rects = {}
  
  for alg in algs:
    alginfo = dsinfo[alg]
    if alg not in series:
      continue
    ymax = max(ymax, max(series[alg]))
    # rects[alg] = axs.errorbar(threads, series[alg], yerr=error[alg],
    rects[alg] = axs.plot(threads, series[alg],
      alpha=opacity,
      color=alginfo.color,
      linewidth=3.0,
      #hatch=hatch[ds],
      linestyle=alginfo.linestyle,
      marker=alginfo.marker,
      markersize=14,
      label=alginfo.name)

  axs.set_ylim(bottom=-0.02*ymax)
  # plt.xticks(threads, threads)
  # axs.set_xlabel('Number of threads')
  # axs.set_ylabel('Throughput (Mop/s)')
  plt.axvline(144, linestyle='--', color='grey') 
  axs.set(xlabel='Number of threads', ylabel='Throughput (Mop/s)')
  legend_x = 1
  legend_y = 0.5 
  # if this_file_name == 'Update_heavy_with_RQ_-_100K_Keys':
  #   plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))



  # plt.legend(framealpha=0.0)
  plt.grid()

  if not paper_ver:
    plt.title(graphtitle)
    plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y))
  
  plt.savefig(outputFile, bbox_inches='tight')

  if paper_ver:
    legend = plt.legend(loc='center left', bbox_to_anchor=(legend_x, legend_y), ncol=8, framealpha=0.0)
    export_legend(legend, 'graphs/' + graph_name + '_legend.pdf')

  
  plt.close('all')

def plot_scalability_graphs(throughput, stddev, threads, ratios, size, algs, graph_name, paper_ver=False):
  for ratio in ratios:
    sufficient_datapoints = False
    for alg in algs:
        num_datapoints = 0
        for th in threads:
            if toString(alg, th, ratio, size) in throughput:
                num_datapoints += 1
        if num_datapoints > 2:
            sufficient_datapoints = True
    if sufficient_datapoints:
        plot_scalability_graph(throughput, stddev, threads, ratio, size, algs, graph_name, paper_ver)


if __name__ == "__main__":
  throughput = {}
  stddev = {}
  threads = []
  ratios = []
  maxkeys = []
  alphas = []
  algs = []

  for filename in input_files:
    readResultsFile(filename, throughput, stddev, threads, ratios, maxkeys, alphas, algs)

  threads.sort()
  ratios.sort()
  maxkeys.sort()
  alphas.sort()

  print('threads: ' + str(threads))
  print('update ratios: ' + str(ratios))
  print('maxkeys: ' + str(maxkeys))
  print('alphas: ' + str(alphas))
  print('algs: ' + str(algs))

  plot_scalability_graph(throughput, stddev, threads, 50, 100000, 0.75, ds_list["try_lock_exp"], "try_vs_strict_lock", True)
  plot_scalability_graph(throughput, stddev, threads, 50, 100000, 0.75, ds_list["trees"], "trees", True)
  plot_scalability_graph(throughput, stddev, threads, 50, 100, 0.75, ds_list["lists"], "lists", True)
  plot_scalability_graph(throughput, stddev, threads, 50, 100000, 0.75, ds_list["rtrees"], "rtrees", True)

  plot_all_graphs(throughput, stddev, threads, ratios, maxkeys, alphas, ds_list["try_lock_exp"], "try_vs_strict_lock")
  plot_all_graphs(throughput, stddev, threads, ratios, maxkeys, alphas, ds_list["trees"], "trees")
  plot_all_graphs(throughput, stddev, threads, ratios, maxkeys, alphas, ds_list["lists"], "lists")
  plot_all_graphs(throughput, stddev, threads, ratios, maxkeys, alphas, ds_list["rtrees"], "rtrees")