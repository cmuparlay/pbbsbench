#include <cstdio>
#include <cstdlib>
#include <cstdint>
// #include <memory>
#include "HNSW.hpp"
#include "benchUtils.h"
using ANN::HNSW;

struct fvec
{
	uint32_t id;
	float *coord;
};

class descr_fvec
{
public:
	typedef fvec type_point;
	static float distance(const type_point &u, const type_point &v, uint32_t dim)
	{
		const auto &uc=u.coord, &vc=v.coord;
		float sum = 0;
		for(uint32_t i=0; i<dim; ++i)
		{
			const auto d = uc[i]-vc[i];
			sum += d*d;
		}
		return sum;
	}
};

int main(int argc, char **argv)
{
	commandLine parameter(argc, argv, "[-q <queryFile>] [-o <outFile>] <inFile>");
	const char* file_in = parameter.getArgument(0);
	const char* file_out = parameter.getOptionValue("-o");
	const char* file_query = parameter.getOptionValue("-q");

	auto to_fvec = [](size_t id, auto begin, auto end){
		typedef typename std::iterator_traits<decltype(begin)>::value_type type_elem;
		if constexpr(std::is_same_v<decltype(begin),ptr_mapped<type_elem,ptr_mapped_src::DISK>>)
		{
			const auto *begin_raw=begin.get(), *end_raw=end.get();
			const auto n = std::distance(begin_raw, end_raw);

			// auto coord = std::make_unique<type_elem[]>(n);
			type_elem *coord = new type_elem[n];
			parlay::parallel_for(0, n, [&](size_t i){
				coord[i] = *(begin_raw+i);
			});

			fvec point;
			point.id = id;
			point.coord = std::move(coord);
			return point;
		}
	};

	parlay::internal::timer t("HNSW", true);
	auto [ps,dim] = parse_vecs<float>(file_in, to_fvec);
	t.next("Read inFile");
	auto [q,_] = parse_vecs<float>(file_query, to_fvec);
	t.next("Read queryFile");

	fputs("Start building HNSW\n", stderr);
	HNSW<descr_fvec> g(dim, ps.begin(), ps.end());
	t.next("Build index");

	// q.size()
	// parlay::parallel_for();
	std::vector<std::vector<fvec*>> res(10);
	parlay::parallel_for(0, 10, [&](size_t  i){
		res[i] = g.search(q[i],10,50);
	});
	t.next("Find neighbors");

	for(uint32_t i=0; i<10; ++i)
	{
		for(uint32_t j=res.size(); j>0; --j)
			printf("%u\t", res[i][j-1]->id);
		putchar('\n');
	}
/*
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
*/
	return 0;
}