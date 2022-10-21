#!/usr/bin/python

bnchmrk="neighbors"
benchmark="Nearest Neighbors"
checkProgram="../bench/neighborsCheck"
dataDir = "/ssd0/ANN/sift10k"

tests = [    
    [1, "siftsmall_base.fvecs", "siftsmall_query.fvecs", "siftsmall_groundtruth.ivecs", 
        "-a 1.2 -R 10 -L 15 -k 10", "-k 10 -r '[1]'"],
    [1, "siftsmall_base.fvecs", "siftsmall_query.fvecs", "siftsmall_groundtruth.ivecs", 
        "-a 1.2 -R 25 -L 35 -k 30", "-k 30 -r '[1, 2, 5, 10, 20]'"],
    [1, "siftsmall_base.fvecs", "siftsmall_query.fvecs", "siftsmall_groundtruth.ivecs", 
        "-a 1.2 -R 50 -L 75 -k 70", "-k 70 -r '[1, 2, 5, 10, 20, 50]'"],
    [1, "siftsmall_base.fvecs", "siftsmall_query.fvecs", "siftsmall_groundtruth.ivecs", 
        "-a 1.2 -R 100 -L 125 -k 100", "-k 100 -r '[1, 2, 5, 10, 20, 50, 75, 100]'"]
    ]


import sys
sys.path.insert(0, 'common')
import runTestsANN
runTestsANN.timeAllArgs(bnchmrk, benchmark, checkProgram, dataDir, tests)

