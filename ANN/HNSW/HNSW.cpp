#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include "HNSW.hpp"
// #include "benchUtils.h"
using namespace std;
using ANN::HNSW;

class fvec
{
public:
	typedef std::vector<double> type_points;
	static double distance(const type_points &u, const type_points &v)
	{
		const auto n = u.size();
		double sum = 0;
		for(uint32_t i=0; i<n; ++i)
		{
			const auto d = u[i]-v[i];
			sum += d*d;
		}
		return sum;
	}
};

int main()
{
	std::vector<std::vector<double>> v;
	v.push_back(std::vector<double>{1,1});
	v.push_back(std::vector<double>{2,3});
	v.push_back(std::vector<double>{4,5});
	HNSW<fvec> g(v.begin(), v.end());
	const auto res = g.search(std::vector<double>{4,4.5}, 10, 50);
	for(const auto &v : res)
	{
		// printf("%f %f\n", v->x, v->y);
		for(const auto c : *v)
			printf("%f ", c);
		putchar('\n');
	}
	return 0;
}