#include "graph.h"
#include "parallel.h"
#include "sequence.h"

// vertexId needs to be signed for the serial version
using vertexId = int;
using edgeId = uint;
using edgeWeight = float;

pbbs::sequence<edgeId> mst(wghEdgeArray<vertexId,edgeWeight> const &E);

