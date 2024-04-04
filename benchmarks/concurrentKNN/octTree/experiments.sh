#!/bin/bash
mkdir -p graphs
#run to quickly test the setup
# python3 run_experiments.py [neighbors_bench,LFKDTree] [36,72,144] [50] [1] 3 [3DinCube_20M]
python3 run_experiments.py [LFKDTree] [144] [0,25,50,75,100] [1] 3 [3DinCube_20M]
python3 run_experiments.py [LFKDTree] [144] [50] [1] 3 [3DinCube_2M,3DinCube_20M,3DinCube_200M,3DinCube_2B]
# python3 run_experiments.py [neighbors_bench] [36,72,144] [50] [1,10] 3 [3DinCube_20M]
# python3 run_experiments.py [working_set_bench,working_set_bench_hoh,working_set_bench_path_copy] [1,2,8,12,36,72,144] [50] [1] 3 [3DinCube_WorkingSet_11M]


# python3 run_experiments.py [neighbors_bench,neighbors_bench_path_copy,neighbors_bench_lockfree,neighbors_bench_hoh] [1,2,4,8,12,36,72,144] [100] [1] 3 [3DinCube_20M]
# python3 run_experiments.py [neighbors_bench,neighbors_bench_path_copy] [1,2,8,12,36,72,144] [100] [1,10] 3 [3DinCube_20M]
# python3 run_experiments.py [neighbors_bench,neighbors_bench_path_copy] [1,2,8,12,36,72,144] [50] [1,10] 3 3DinCube_20M
# python3 run_experiments.py [neighbors_bench,neighbors_bench_path_copy] [1,2,8,12,36,72,144] [10] [1,10] 3 3DinCube_20M



# python3 run_experiments.py [neighbors_bench,neighbors_bench_path_copy] [1,2,4,8,12,36,72,144] [50] [1,10] 2 [2DinCube_20M]



# python3 run_experiments.py [neighbors_bench,neighbors_bench_path_copy] [144] [0,25,50,75,100] [1] 3 [3Dplummer_20M]
# python3 run_experiments.py [neighbors_bench,neighbors_bench_path_copy] [144] [0,25,50,75,100] [1] 3 [3DinCube_20M]

# python3 run_experiments.py [neighbors_bench,neighbors_bench_path_copy] [1,2,8,12,36,72,144] [10] [1,10] 3 [3Dplummer_20M]

# python3 run_experiments.py [neighbors_bench,neighbors_bench_path_copy] [1,2,8,12,36,72,144] [50] [1,10] 3 [lucy3D_14M]

# python3 run_experiments.py [neighbors_bench,neighbors_bench_path_copy] [1,2,8,12,36,72,144] [50] [1,10] 3 [thai_statue5M]

# python3 run_experiments.py [range_bench,range_bench_path_copy] [1,2,8,12,36,72,144] [10] [.014,.0176] 3 [3Dplummer_20M]

# python3 run_experiments.py [neighbors_bench] [144] [50] [1] 3 [3DinCube_2M,3DinCube_20M,3DinCube_200M,3DinCube_2B]

##OVERSUBSCRIPTION

# python3 run_experiments.py [neighbors_bench,neighbors_bench_lockfree,neighbors_bench_path_copy,neighbors_bench_path_copy_lockfree] [1,2,8,12,36,72,144,196,288,432,500] [50] [1] 3 [3DinCube_20M]
# python3 run_experiments.py [working_set_bench,working_set_bench_lockfree] [1,2,8,12,36,72,144,196,288,432,500] [50] [1] 3 [3DinCube_WorkingSet_11M]




