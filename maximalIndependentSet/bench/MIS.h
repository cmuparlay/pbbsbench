#include "graph.h"

using vertexId = uint;
using edgeId = uint;
using Graph = graph<vertexId,edgeId>;

pbbs::sequence<char> maximalIndependentSet(Graph G);

