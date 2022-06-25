#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <string>
#include <map>
// #include <memory>
#include <H5Cpp.h>
#include "HNSW.hpp"
#include "dist_ang.hpp"
using ANN::HNSW;

// it will change the file string
auto load_fvec(char *file)
{
	char *spec_input = std::strchr(file, ':');
	if(spec_input==nullptr)
	{
		fputs("Unrecognized file spec",stderr);
		return std::make_pair(parlay::sequence<fvec>(),uint32_t(0));
	}
	
	*(spec_input++) = '\0';
	parlay::sequence<fvec> ps;
	uint32_t dim = 0;
	if(spec_input[0]=='/')
		std::tie(ps,dim) = ps_from_HDF5(file, spec_input);
	else if(std::strcmp(spec_input,"fvec"))
		std::tie(ps,dim) = ps_from_SIFT(file);
	else fputs("Unsupported file spec",stderr);

	return std::make_pair(ps,dim);
}

// it will change the file string
auto load_ivec(char *file)
{
	char *spec_input = std::strchr(file, ':');
	if(spec_input==nullptr)
	{
		fputs("Unrecognized file spec",stderr);
		return std::make_pair(parlay::sequence<uint32_t*>(),uint32_t(0));
	}
	
	*(spec_input++) = '\0';
	parlay::sequence<uint32_t*> vec;
	uint32_t bound1 = 0;
	if(spec_input[0]=='/')
	{
		auto [buffer_ptr,bound] = read_array_from_HDF5<uint32_t>(file, spec_input);
		bound1 = bound[1];
		auto *buffer = buffer_ptr.release();

		vec.resize(bound[0]);
		parlay::parallel_for(0, bound[0], [&](uint32_t i){
			vec[i] = &buffer[i*bound1];
		});
	}
	else if(std::strcmp(spec_input,"fvec"))
		std::tie(vec,bound1) = parse_vecs<uint32_t>(file, [](size_t, auto begin, auto end){
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
	else fputs("Unsupported file spec",stderr);

	return std::make_pair(vec,bound1);
}

typedef void (*type_func)(HNSW<descr_fvec> &, commandLine, parlay::internal::timer&);

void output_deg(HNSW<descr_fvec> &g, commandLine param, parlay::internal::timer &t)
{
	if(param.getOption("-?"))
	{
		printf(__func__);
		puts(": -f <outFile> [-l <level>=0]");
		return;
	};
	const char* outfile = param.getOptionValue("-f");
	const uint32_t level = param.getOptionIntValue("-l", 0);
	auto deg = g.get_deg(level);
	FILE *file_deg = fopen(outfile, "w");
	if(!file_deg)
	{
		fputs("fail to create the file\n", stderr);
		exit(2);
	}
	for(auto e : deg)
	{
		fprintf(file_deg, "%u\n", e);
	}
	fclose(file_deg);
}

void output_recall(HNSW<descr_fvec> &g, commandLine param, parlay::internal::timer &t)
{
	if(param.getOption("-?"))
	{
		printf(__func__);
		puts(
			"[-q <queryFile>] [-g <groundtruthFile>]"
			"-ef <ef_query> [-r <recall@R>=1] [-k <numQuery>=all]"
		);
		return;
	};
	char* file_query = param.getOptionValue("-q");
	char* file_groundtruth = param.getOptionValue("-g");
	auto [q,_] = load_fvec(file_query);
	t.next("Read queryFile");

	uint32_t cnt_rank_cmp = param.getOptionIntValue("-r", 1);
	const uint32_t ef = param.getOptionIntValue("-ef", cnt_rank_cmp*50);
	const uint32_t cnt_pts_query = param.getOptionIntValue("-k", q.size());

	std::vector<std::vector<std::pair<uint32_t,double>>> res(cnt_pts_query);
	parlay::parallel_for(0, cnt_pts_query, [&](size_t i){
		res[i] = g.search(q[i], cnt_rank_cmp, ef);
	});
	t.next("Find neighbors");

	auto [gt,rank_max] = load_ivec(file_groundtruth);

	if(rank_max<cnt_rank_cmp)
		cnt_rank_cmp = rank_max;
	uint32_t cnt_all_shot = 0;
	printf("measure recall@%u\n", cnt_rank_cmp);
	for(uint32_t i=0; i<cnt_pts_query; ++i)
	{
		uint32_t cnt_shot = 0;
		for(uint32_t j=0; j<cnt_rank_cmp; ++j)
			if(std::find_if(res[i].begin(),res[i].end(),[&](const std::pair<uint32_t,double> &p){
				return p.first==gt[i][j];}) != res[i].end())
			{
				cnt_shot++;
			}
		printf("#%u:\t%u (%.2f)[%lu]", i, cnt_shot, float(cnt_shot)/cnt_rank_cmp, res[i].size());
		if(cnt_shot==cnt_rank_cmp)
		{
			cnt_all_shot++;
		}
		putchar('\n');
	}
	printf("#all shot: %u (%.2f)\n", cnt_all_shot, float(cnt_all_shot)/cnt_pts_query);
}

void output_neighbor(HNSW<descr_fvec> &g, commandLine param, parlay::internal::timer &t)
{
	if(param.getOption("-?"))
	{
		printf(__func__);
		puts("[-q <queryFile>] [<ef_query> <recall> <begin> <end> <stripe>]...");
		return;
	};
	char* file_query = param.getOptionValue("-q");
	auto [q,_] = load_fvec(file_query);

	puts("Please input <ef_query> <recall> <begin> <end> <stripe> in order");
	while(true)
	{
		uint32_t ef, recall, begin, end, stripe;
		scanf("%u%u%u%u%u", &ef, &recall, &begin, &end, &stripe);
		for(uint32_t i=begin; i<end; i+=stripe)
		{
			auto res = g.search_ex(q[i], recall, ef);
			printf("Neighbors of %u\n", i);
			for(auto it=res.crbegin(); it!=res.crend(); ++it)
			{
				const auto [id,dep,dis] = *it;
				printf("  [%u]\t%u\t%.6f\n", dep, id, dis);
			}
			putchar('\n');
		}
	}
}

void count_number(HNSW<descr_fvec> &g, commandLine param, parlay::internal::timer &t)
{
	if(param.getOption("-?"))
	{
		printf(__func__);
		return;
	};
	std::map<uint32_t,uint32_t> cnt;
	for(const auto *p : g.node_pool)
		cnt[p->level]++;
	
	uint32_t sum = 0;
	for(int i=cnt.rbegin()->first; i>=0; --i)
		printf("#nodes in lev. %d: %u (%u)\n", i, sum+=cnt[i], cnt[i]);
}

int main(int argc, char **argv)
{
	for(int i=0; i<argc; ++i)
		printf("%s ", argv[i]);
	putchar('\n');

	commandLine param(argc, argv, 
		"-in <inFile> -mod <modelFile> -func <function>"
	);
	char* file_in = param.getOptionValue("-in");
	const char* file_model = param.getOptionValue("-mod");
	const char* func = param.getOptionValue("-func");

	parlay::internal::timer t("HNSW", true);
	auto [ps,_] = load_fvec(file_in);
	t.next("Read inFile");

	fputs("Start building HNSW\n", stderr);
	HNSW<descr_fvec> g(file_model, [&](uint32_t i){
		return ps[i];
	});
	t.next("Build index");

	std::map<std::string,type_func> list_func;
	list_func["deg"] = output_deg;
	list_func["recall"] = output_recall;
	list_func["neighbor"] = output_neighbor;
	list_func["count"] = count_number;
	auto it_func = list_func.find(func);
	if(it_func!=list_func.end())
		it_func->second(g, param, t);
	else
		fprintf(stderr, "Cannot find function '%s'\n", func);

	return 0;
}
