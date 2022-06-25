#ifndef __H5_OPS_HPP__
#define __H5_OPS_HPP__

#include <cstdio>
#include <array>
#include <memory>
#include <type_traits>
#include <H5Cpp.h>

// read a 2D array from H5 file and return a 1D array
template<typename T>
std::pair<std::unique_ptr<T[]>,std::array<uint32_t,2>> read_array_from_HDF5(const char *file, const char *dir)
{
	H5::H5File file_h5(file, H5F_ACC_RDONLY);
	H5::DataSet dset = file_h5.openDataSet(dir);
	H5::DataSpace dspace = dset.getSpace();
	hsize_t bound[2];
	dspace.getSimpleExtentDims(bound);
	fprintf(stderr, "%s: [%llu,%llu]\n", dir, bound[0], bound[1]);

	hsize_t bound_1d = bound[0]*bound[1];
	auto buffer = std::make_unique<T[]>(bound_1d);

	static_assert(std::is_arithmetic<T>::value, "Unsupported type");
	auto type_H5data = std::is_integral_v<T>?
		H5::PredType::NATIVE_UINT32: H5::PredType::NATIVE_FLOAT;
	dset.read(buffer.get(), type_H5data, H5::DataSpace(1,&bound_1d,NULL), dspace);

	return {std::move(buffer), std::array<uint32_t,2>{{bound[0],bound[1]}}};
}

#endif // __H5_OPS_HPP__
