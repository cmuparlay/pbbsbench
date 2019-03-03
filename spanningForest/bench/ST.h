#include "graph.h"
#include "sequence.h"

// vertexId needs to be signed for the serial version
using vertexId = int;
using edgeId = uint;

pbbs::sequence<edgeId> st(edgeArray<vertexId> const &EA);
