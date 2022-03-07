#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <algorithm>
// #include "err.h"
#include "benchUtils.h"

template<typename T>
struct vec
{
	int id;
	parlay::slice<T*, T*> coordinates;

	vec() :
		coordinates(parlay::make_slice<T*,T*>(nullptr,nullptr)){
	}
};

int main(int argc, char **argv)
{
	if(argc<2)
	{
		fprintf(stderr, "Usage: %s <fvecFile> [cntPointPrinted] [cntDimPrinted]\n", argv[0]);
		exit(1);
	}

	const auto points = parse_vecs<float>(argv[1], [](size_t id, auto begin, auto end){
		typedef typename std::iterator_traits<decltype(begin)>::value_type type_elem;
		static_assert(std::is_same_v<decltype(begin),ptr_mapped<type_elem,ptr_mapped_src::DISK>>);

		vec<type_elem> point;
		point.id = id;
		point.coordinates = parlay::make_slice(begin.get(), end.get());
		return point;
	});

	uint32_t cnt_ps=points.size(), cnt_dim=points[0].coordinates.size();
	if(argc>=3)
	{
		const uint32_t cnt_ps_cmd = atoi(argv[2]);
		if(cnt_ps_cmd<cnt_ps)
			cnt_ps = cnt_ps_cmd;
	}
	if(argc>=4)
	{
		const uint32_t cnt_dim_cmd = atoi(argv[3]);
		if(cnt_dim_cmd<cnt_dim)
			cnt_dim = cnt_dim_cmd;
	}

	for(int i=0; i<cnt_ps; ++i)
	{
		const auto &p = points[i];
		printf("[#%d]\t", p.id);
		for(int j=0; j<cnt_dim; ++j)
			printf(" %8.2f", p.coordinates[j]);
		putchar('\n');
	}
	return 0;
}