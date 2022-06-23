#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <string>
#include <map>
// #include <memory>
#include <H5Cpp.h>
#include "HNSW.hpp"
#include "benchUtils.h"
using ANN::HNSW;

struct fvec
{
	uint32_t id;
	float *coord;
};

class descr_fvec_l2
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

class descr_fvec
{
public:
	typedef fvec type_point;
	static double distance(const type_point &u, const type_point &v, uint32_t dim)
	{
		const auto &uc=u.coord, &vc=v.coord;
		double dot=0, nu=0, nv=0;
		for(uint32_t i=0; i<dim; ++i)
		{
			nu += double(uc[i])*uc[i];
			nv += double(vc[i])*vc[i];
			dot += double(uc[i])*vc[i];
		}
		return 1-dot/(sqrt(nu)*sqrt(nv));
	}

	static auto get_id(const type_point &u)
	{
		return u.id;
	}
};

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
	const char* file_query = param.getOptionValue("-in");
	/*
	const char* file_query = param.getOptionValue("-q");
	const char* file_groundtruth = param.getOptionValue("-g");
	auto [q,_] = parse_vecs<float>(file_query, to_fvec);
	t.next("Read queryFile");
	*/
	H5::H5File file(file_query, H5F_ACC_RDONLY);
	H5::DataSet dset_test = file.openDataSet("/test");
	H5::DataSpace dspace_test = dset_test.getSpace();
	hsize_t bound[2];
	dspace_test.getSimpleExtentDims(bound);
	printf("/test: [%llu,%llu]\n", bound[0], bound[1]);

	hsize_t bound_1d = bound[0]*bound[1];
	auto buf_test = new float[bound_1d];
	dset_test.read(buf_test, H5::PredType::NATIVE_FLOAT, H5::DataSpace(1,&bound_1d,NULL), dspace_test);
	parlay::sequence<fvec> q(bound[0]);
	parlay::parallel_for(0, bound[0], [&](uint32_t i){
		q[i] = fvec{i, &buf_test[i*bound[1]]};
	});
	dspace_test.close();
	dset_test.close();

	const char* ef = param.getOptionValue("-ef");
	uint32_t cnt_rank_cmp = param.getOptionIntValue("-r", 1);
	const uint32_t cnt_pts_query = param.getOptionIntValue("-k", q.size());

	std::vector<std::vector<std::pair<uint32_t,double>>> res(cnt_pts_query);
	parlay::parallel_for(0, cnt_pts_query, [&](size_t i){
		res[i] = g.search(q[i], cnt_rank_cmp, atoi(ef));
	});
	t.next("Find neighbors");
/*
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
*/
	H5::DataSet dset_gt = file.openDataSet("/neighbors");
	H5::DataSpace dspace_gt = dset_gt.getSpace();
	dspace_gt.getSimpleExtentDims(bound);
	auto rank_max = bound[1];

	bound_1d = bound[0]*bound[1];
	auto buf_gt = new uint32_t[bound_1d];
	dset_gt.read(buf_gt, H5::PredType::NATIVE_UINT32, H5::DataSpace(1,&bound_1d,NULL), dspace_gt);
	dspace_gt.close();
	dset_gt.close();
	file.close();

	if(rank_max<cnt_rank_cmp)
		cnt_rank_cmp = rank_max;
	uint32_t cnt_all_shot = 0;
	printf("measure recall@%u\n", cnt_rank_cmp);
	for(uint32_t i=0; i<cnt_pts_query; ++i)
	{
		uint32_t cnt_shot = 0;
		for(uint32_t j=0; j<cnt_rank_cmp; ++j)
			if(std::find_if(res[i].begin(),res[i].end(),[&](const std::pair<uint32_t,double> &p){
			//	return p.first==gt[i][j];}) != res[i].end())
				return p.first==buf_gt[i*rank_max+j];}) != res[i].end())
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
//	const char* file_query = param.getOptionValue("-q");
	const char* file_query = param.getOptionValue("-in");
//	auto [q,_] = parse_vecs<float>(file_query, to_fvec);

	H5::H5File file(file_query, H5F_ACC_RDONLY);
	H5::DataSet dset_test = file.openDataSet("/test");
	H5::DataSpace dspace_test = dset_test.getSpace();
	hsize_t bound[2];
	dspace_test.getSimpleExtentDims(bound);
	printf("/test: [%llu,%llu]\n", bound[0], bound[1]);

	hsize_t bound_1d = bound[0]*bound[1];
	auto buf_test = new float[bound_1d];
	dset_test.read(buf_test, H5::PredType::NATIVE_FLOAT, H5::DataSpace(1,&bound_1d,NULL), dspace_test);
	parlay::sequence<fvec> q(bound[0]);
	parlay::parallel_for(0, bound[0], [&](uint32_t i){
		q[i] = fvec{i, &buf_test[i*bound[1]]};
	});
	dspace_test.close();
	dset_test.close();
	file.close();

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
	const char* file_in = param.getOptionValue("-in");
	const char* file_model = param.getOptionValue("-mod");
	const char* func = param.getOptionValue("-func");

	parlay::internal::timer t("HNSW", true);
//	auto [ps,dim] = parse_vecs<float>(file_in, to_fvec);

	H5::H5File file(file_in, H5F_ACC_RDONLY);
	H5::DataSet dset_train = file.openDataSet("/train");
	H5::DataSpace dspace_train = dset_train.getSpace();
	hsize_t bound[2];
	dspace_train.getSimpleExtentDims(bound);
	const auto dim = bound[1];
	printf("/train: [%llu,%llu]\n", bound[0], bound[1]);

	hsize_t bound_1d = bound[0]*bound[1];
	auto buf_train = new float[bound_1d];
	dset_train.read(buf_train, H5::PredType::NATIVE_FLOAT, H5::DataSpace(1,&bound_1d,NULL), dspace_train);
	parlay::sequence<fvec> ps(bound[0]);
	parlay::parallel_for(0, bound[0], [&](uint32_t i){
		ps[i] = fvec{i, &buf_train[i*dim]};
	});
	dspace_train.close();
	dset_train.close();
	file.close();
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
