#include "common/graph.h"
#include "parlay/primitives.h"

// vertexId needs to be signed for the serial version
using vertexId = int;
using edgeId = uint;
using edgeWeight = float;

parlay::sequence<edgeId> mst(wghEdgeArray<vertexId,edgeWeight> &E);


