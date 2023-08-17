

# lock-free, no path copying
g++ -DHOMEGROWN -pthread -mcx16 -O0 -g -std=c++17  -I . -DVersioned -DLazyStamp -include neighbors.h -o neighbors_lf ../bench/neighborsTime.C -DHOMEGROWN -pthread -ldl -L/usr/local/lib -ljemalloc

# lock-based, no path copying
g++ -DHOMEGROWN -pthread -mcx16 -O0 -g -std=c++17  -I . -DNoHelp -DVersioned -DLazyStamp -include neighbors.h -o neighbors ../bench/neighborsTime.C -DHOMEGROWN -pthread -ldl -L/usr/local/lib -ljemalloc

# lock-based, no path copying
g++ -DHOMEGROWN -pthread -mcx16 -O0 -g -std=c++17  -I . -DPathCopy -DNoHelp -DVersioned -DLazyStamp -include neighbors.h -o neighbors_path_copy ../bench/neighborsTime.C -DHOMEGROWN -pthread -ldl -L/usr/local/lib -ljemalloc

# lock-free, path copying
g++ -DHOMEGROWN -pthread -mcx16 -O0 -g -std=c++17  -I . -DPathCopy -DVersioned -DLazyStamp -include neighbors.h -o neighbors_path_copy_lf ../bench/neighborsTime.C -DHOMEGROWN -pthread -ldl -L/usr/local/lib -ljemalloc