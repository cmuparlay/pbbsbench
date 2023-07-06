#!/bin/bash

set -o xtrace

# lock-based
g++ -DHOMEGROWN -pthread -mcx16 -O3 -std=c++17 -DNDEBUG -I . -DNoHelp -DVersioned -DLazyStamp   -include range_bench.h -o range_bench ../bench/rangeTime.C -DHOMEGROWN -pthread -ldl -L/usr/local/lib -ljemalloc


# lock-based, path copying

g++ -DHOMEGROWN -pthread -mcx16 -O3 -std=c++17 -DNDEBUG -I . -DNoHelp -DVersioned -DPathCopy -DLazyStamp -include range_bench.h -o range_bench_path_copy ../bench/rangeTime.C -DHOMEGROWN -pthread -ldl -L/usr/local/lib -ljemalloc
