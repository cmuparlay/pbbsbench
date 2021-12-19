#include <cstdio>
#include <cstdlib>
#include "HNSW.hpp"
// #include "benchUtils.h"
using namespace std;
using ANN::HNSW;

int main()
{
	std::vector<ANN::point> v;
	v.push_back({1,1});
	v.push_back({2,3});
	v.push_back({4,5});
	HNSW<ANN::point> g(v.begin(), v.end());
	const auto res = g.search(ANN::point{4,4.5}, 10, 50);
	for(const auto &v : res)
		printf("%f %f\n", v->x, v->y);
	return 0;
}