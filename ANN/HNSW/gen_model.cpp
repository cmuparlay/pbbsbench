#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <algorithm>
// #include <memory>
#include <tuple>
#include "HNSW.hpp"
// #include "benchUtils.h" // fix multiple inclusion
// #include "dist_l2.hpp"
#include "dist_ang.hpp"
using ANN::HNSW;

int main(int argc, char **argv)
{
	for(int i=0; i<argc; ++i)
		printf("%s ", argv[i]);
	putchar('\n');

	commandLine parameter(argc, argv, 
		"-n <numInput> -ml <m_l> -m <m> "
		"-efc <ef_construction> -alpha <alpha> -r <recall@R> [-b <batchBase>]"
		"-in <inFile> -out <modelFile>"
	);
	char* file_in = parameter.getOptionValue("-in");
	const char* file_out = parameter.getOptionValue("-out");
	const uint32_t cnt_points = parameter.getOptionLongValue("-n", 0);
	const float m_l = parameter.getOptionDoubleValue("-ml", 0.36);
	const uint32_t m = parameter.getOptionIntValue("-m", 40);
	const uint32_t efc = parameter.getOptionIntValue("-efc", 60);
	const float alpha = parameter.getOptionDoubleValue("-alpha", 1);
	const float batch_base = parameter.getOptionDoubleValue("-b", 2);
	const bool do_fixing = !!parameter.getOptionIntValue("-f", 0);

	if(file_in==nullptr || file_out==nullptr)
		return fputs("in/out files are not indicated\n",stderr), 1;
	
	char *spec_input = std::strchr(file_in, ':');
	if(spec_input==nullptr)
		return fputs("Unrecognized file spec",stderr), 2;
	
	parlay::internal::timer t("HNSW", true);

	*(spec_input++) = '\0';
	parlay::sequence<fvec> ps;
	uint32_t dim;
	if(spec_input[0]=='/')
		std::tie(ps,dim) = ps_from_HDF5(file_in, spec_input);
	else if(std::strcmp(spec_input,"fvec"))
		std::tie(ps,dim) = ps_from_SIFT(file_in);
	else return fputs("Unsupported file spec",stderr), 2;

	t.next("Read inFile");

	fputs("Start building HNSW\n", stderr);
	HNSW<descr_fvec> g(
		ps.begin(), ps.begin()+cnt_points, dim,
		m_l, m, efc, alpha, batch_base, do_fixing
	);
	t.next("Build index");

	g.save(file_out);

	t.next("Write to the file");
	return 0;
}
