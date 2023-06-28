#!/bin/bash

set -o xtrace

# lock-based
g++ -DHOMEGROWN -pthread -mcx16 -O3 -std=c++17 -DNDEBUG -I . -DNoHelp -DVersioned -DLazyStamp   -include neighbors_bench.h -o neighbors_bench ../bench/neighborsTime.C -DHOMEGROWN -pthread -ldl -L/usr/local/lib -ljemalloc

# lock-free (still need to fix something)
# g++ -DHOMEGROWN -pthread -mcx16 -O3 -std=c++17 -DNDEBUG -I . -DVersioned -DLazyStamp   -include neighbors_bench.h -o neighbors_bench ../bench/neighborsTime.C -DHOMEGROWN -pthread -ldl -L/usr/local/lib -ljemalloc