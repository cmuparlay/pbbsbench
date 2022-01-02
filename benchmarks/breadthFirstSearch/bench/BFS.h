#include "common/graph.h"
#include "parlay/sequence.h"

// vertexId needs to be signed
using vertexId = int;
using edgeId = uint;
using Graph = graph<vertexId,edgeId>;

// returns a parent sequence where -1 means it was not visited,
// and the start points to itself.
parlay::sequence<vertexId> BFS(vertexId start, const Graph &G, bool verbose);
