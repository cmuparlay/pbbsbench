#ifndef __DIST_L2_HPP__
#define __DIST_L2_HPP__

#include "type_point.hpp"

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


#endif // _DIST_L2_HPP_
