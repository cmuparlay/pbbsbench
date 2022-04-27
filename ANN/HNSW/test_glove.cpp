#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <algorithm>
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

int main(int argc, char **argv)
{
	for(int i=0; i<argc; ++i)
		printf("%s ", argv[i]);
	putchar('\n');

	commandLine parameter(argc, argv, 
		"-n <numInput> -k <numQuery> -ml <m_l> -m <m> "
		"-efc <ef_construction> -alpha <alpha> -ef <ef_query> -r <recall@R> [-b <batchBase>]"
		"<inFile>"
	);
	const char* file_in = parameter.getArgument(0);
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

	parlay::internal::timer t("HNSW", true);

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
	t.next("Read inFile");

	fputs("Start building HNSW\n", stderr);
	HNSW<descr_fvec> g(
		ps.begin(), ps.begin()+atoi(cnt_pts_input), dim,
		atof(m_l), atoi(m), atoi(efc), atof(alpha), atof(batch_base), !!atoi(do_fixing)
	);
	t.next("Build index");

	H5::DataSet dset_test = file.openDataSet("/test");
	H5::DataSpace dspace_test = dset_test.getSpace();
	dspace_test.getSimpleExtentDims(bound);

	bound_1d = bound[0]*bound[1];
	auto buf_test = new float[bound_1d];
	dset_test.read(buf_test, H5::PredType::NATIVE_FLOAT, H5::DataSpace(1,&bound_1d,NULL), dspace_test);
	parlay::sequence<fvec> q(bound[0]);
	parlay::parallel_for(0, bound[0], [&](uint32_t i){
		q[i] = fvec{i, &buf_test[i*bound[1]]};
	});
	dspace_test.close();
	dset_test.close();
	t.next("Read queryFile");

	const uint32_t cnt_pts_query_val = atoi(cnt_pts_query);
	uint32_t cnt_rank_cmp_val = atoi(cnt_rank_cmp);
	std::vector<std::vector<fvec*>> res(cnt_pts_query_val);
	parlay::parallel_for(0, cnt_pts_query_val, [&](size_t i){
		res[i] = g.search(q[i], cnt_rank_cmp_val, atoi(ef));
	});
	t.next("Find neighbors");

	H5::DataSet dset_gt = file.openDataSet("/neighbors");
	H5::DataSpace dspace_gt = dset_gt.getSpace();
	dspace_gt.getSimpleExtentDims(bound);
	auto rank_max = bound[1];

	bound_1d = bound[0]*bound[1];
	auto buf_gt = new uint32_t[bound_1d];
	dset_gt.read(buf_gt, H5::PredType::NATIVE_UINT32, H5::DataSpace(1,&bound_1d,NULL), dspace_gt);
	dspace_gt.close();
	dset_gt.close();

	if(rank_max<cnt_rank_cmp_val)
		cnt_rank_cmp_val = rank_max;
	uint32_t cnt_all_shot = 0;
	printf("measure recall@%u\n", cnt_rank_cmp_val);
	for(uint32_t i=0; i<cnt_pts_query_val; ++i)
	{
		uint32_t cnt_shot = 0;
		for(uint32_t j=0; j<cnt_rank_cmp_val; ++j)
			if(std::find_if(res[i].begin(),res[i].end(),[&](fvec *p){
				return p->id==buf_gt[i*rank_max+j];}) != res[i].end())
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
			printf(" %u[%.4f]", r->id, descr_fvec::distance(q[i],ps[r->id],dim));
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
	file.close();
	return 0;
}
