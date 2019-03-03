#include "graph.h"
#include "parallel.h"
#include "sequence.h"

using vertexId = uint;
using edgeId = uint;
using edgeWeight = float;

pbbs::sequence<vertexId> mst(wghEdgeArray<vertexId,edgeWeight> const &E);

