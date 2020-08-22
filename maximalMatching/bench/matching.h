#include "common/graph.h"

using vertexId = uint;
using edgeId = uint;
using edges = edgeArray<vertexId>;

parlay::sequence<edgeId> maximalMatching(edges const &EA);

