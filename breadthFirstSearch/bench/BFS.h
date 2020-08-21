#include "common/graph.h"

// vertexId needs to be signed
using vertexId = int;
using edgeId = uint;
using Graph = graph<vertexId,edgeId>;

std::pair<vertexId,size_t> BFS(vertexId start, Graph &G);
