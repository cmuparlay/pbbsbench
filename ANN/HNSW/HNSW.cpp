#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
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

	static auto get_id(const type_point &u)
	{
		return u.id;
	}
};

int main(int argc, char **argv)
{
	for(int i=0; i<argc; ++i)
		printf("%s ", argv[i]);
	putchar('\n');

	commandLine parameter(argc, argv, 
		"-n <numInput> -k <numQuery> -ml <m_l> -m <m> "
		"-efc <ef_construction> -alpha <alpha> -ef <ef_query> -r <recall@R> [-b <batchBase>]"
		"[-q <queryFile>] [-g <groundtruthFile>] <inFile>"
	);
	const char* file_in = parameter.getArgument(0);
	const char* file_query = parameter.getOptionValue("-q");
	const char* file_groundtruth = parameter.getOptionValue("-g");
	const char* cnt_pts_input = parameter.getOptionValue("-n");
	const char* cnt_pts_query = parameter.getOptionValue("-k");
	const char* m_l = parameter.getOptionValue("-ml");
	const char* m = parameter.getOptionValue("-m");
	const char* efc = parameter.getOptionValue("-efc");
	const char* alpha = parameter.getOptionValue("-alpha");
	const char* ef = parameter.getOptionValue("-ef");
	const char* cnt_rank_cmp = parameter.getOptionValue("-r");
	const char* batch_base = parameter.getOptionValue("-b");
	const char* do_fixing = parameter.getOptionValue("-f");

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

	fputs("Start building HNSW\n", stderr);
	HNSW<descr_fvec> g("model.bin", [&](uint32_t i){
		return ps[i];
	});
	/*
	HNSW<descr_fvec> g(
		ps.begin(), ps.begin()+atoi(cnt_pts_input), dim,
		atof(m_l), atoi(m), atoi(efc), atof(alpha), atof(batch_base), !!atoi(do_fixing)
	);
	*/
	t.next("Build index");
	// g.save("model.bin");

	if(!file_query) return 0;

	auto [q,_] = parse_vecs<float>(file_query, to_fvec);
	t.next("Read queryFile");

	const uint32_t cnt_pts_query_val = atoi(cnt_pts_query);
	uint32_t cnt_rank_cmp_val = atoi(cnt_rank_cmp);
	std::vector<std::vector<std::pair<uint32_t,double>>> res(cnt_pts_query_val);
	parlay::parallel_for(0, cnt_pts_query_val, [&](size_t i){
		res[i] = g.search(q[i], cnt_rank_cmp_val, atoi(ef));
	});
	t.next("Find neighbors");

	if(!file_groundtruth) return 0;

	auto [gt,rank_max] = parse_vecs<uint32_t>(file_groundtruth, [](size_t, auto begin, auto end){
		typedef typename std::iterator_traits<decltype(begin)>::value_type type_elem;
		if constexpr(std::is_same_v<decltype(begin),ptr_mapped<type_elem,ptr_mapped_src::DISK>>)
		{
			const auto *begin_raw=begin.get(), *end_raw=end.get();
			const auto n = std::distance(begin_raw, end_raw);

			type_elem *id = new type_elem[n];
			parlay::parallel_for(0, n, [&](size_t i){
				id[i] = static_cast<type_elem>(*(begin_raw+i));
			});
			return id;
		}
	});

	if(rank_max<cnt_rank_cmp_val)
		cnt_rank_cmp_val = rank_max;
	uint32_t cnt_all_shot = 0;
	printf("measure recall@%u\n", cnt_rank_cmp_val);
	for(uint32_t i=0; i<cnt_pts_query_val; ++i)
	{
		uint32_t cnt_shot = 0;
		for(uint32_t j=0; j<cnt_rank_cmp_val; ++j)
			if(std::find_if(res[i].begin(),res[i].end(),[&](const std::pair<uint32_t,double> &p){
				return p.first==gt[i][j];}) != res[i].end())
			{
				cnt_shot++;
			}
		//if(cnt_rank_cmp_val==1&&fabs(descr_fvec::distance(q[i],ps[gt[i][0]],dim)-descr_fvec::distance(q[i],ps[res[i][0]->id],dim))<1e-6) cnt_shot=1;
		printf("#%u:\t%u (%.2f)[%lu]", i, cnt_shot, float(cnt_shot)/cnt_rank_cmp_val, res[i].size());
		if(cnt_shot==cnt_rank_cmp_val)
		{
			cnt_all_shot++;
		}
		/*
		for(const auto *r : res[i])
		{
			printf(" %u", r->id);
		}
		*/
		putchar('\n');
	}
	printf("#all shot: %u (%.2f)\n", cnt_all_shot, float(cnt_all_shot)/cnt_pts_query_val);
/*
	for(uint32_t i=0; i<10; ++i)
	{
		for(uint32_t j=res.size(); j>0; --j)
			printf("%u\t", res[i][j-1]->id);
		putchar('\n');
	}
*/
	return 0;
}
