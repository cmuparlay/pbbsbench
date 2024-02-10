#!/bin/bash

set -o xtrace

# lock-based
g++ -DHOMEGROWN -pthread -mcx16 -O3 -std=c++17 -DNDEBUG -I . -DNoHelp -DVersioned -DHWStamp -include neighbors_bench.h -o neighbors_bench ../bench/neighborsTime.C -DHOMEGROWN -pthread -ldl -L/usr/local/lib -ljemalloc

#lock-based, hand over hand
g++ -DHOMEGROWN -pthread -mcx16 -O3 -std=c++17 -DNDEBUG -I . -DNoHelp -DVersioned -DHWStamp -DHandOverHand -include neighbors_bench.h -o neighbors_bench_hoh ../bench/neighborsTime.C -DHOMEGROWN -pthread -ldl -L/usr/local/lib -ljemalloc

# lock-based, path copying
g++ -DHOMEGROWN -pthread -mcx16 -O3 -std=c++17 -DNDEBUG -I . -DNoHelp -DVersioned -DPathCopy -DHWStamp -include neighbors_bench.h -o neighbors_bench_path_copy ../bench/neighborsTime.C -DHOMEGROWN -pthread -ldl -L/usr/local/lib -ljemalloc

# lock-free
g++ -DHOMEGROWN -pthread -mcx16 -O3 -std=c++17 -DNDEBUG -I . -DVersioned -DLazyStamp   -include neighbors_bench.h -o neighbors_bench_lockfree ../bench/neighborsTime.C -DHOMEGROWN -pthread -ldl -L/usr/local/lib -ljemalloc