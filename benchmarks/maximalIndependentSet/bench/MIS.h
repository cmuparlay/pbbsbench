#include "common/graph.h"

using vertexId = uint;
using edgeId = uint;
using Graph = graph<vertexId,edgeId>;

parlay::sequence<char> maximalIndependentSet(Graph const &G);

