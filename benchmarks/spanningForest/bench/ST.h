#include "common/graph.h"

// vertexId needs to be signed for the serial version
using vertexId = int;
using edgeId = uint;

parlay::sequence<edgeId> st(edgeArray<vertexId> const &EA);
