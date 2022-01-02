// Boost.Geometry Index
//
// Modified from Boost Quickbook Examples
//
// Copyright (c) 2011-2013 Adam Wulkiewicz, Lodz, Poland.
//
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>

#include <boost/geometry/index/rtree.hpp>

// to store queries results
#include <vector>

// just for output
#include <iostream>
#include <omp.h>
#include <boost/foreach.hpp>
#include "parlay/primitives.h"
#include "parlay/internal/get_time.h"
#include "range.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
using namespace std;

typedef coord data_type;

typedef pair<data_type, data_type> point_type;


long range(Points const &pts, Queries const &queries, bool verbose) 
{
	
	std::cout << "start running ..." << std::endl;
    typedef bg::model::point<data_type, 2, bg::cs::cartesian> point;
    typedef bg::model::segment<point> segment;
    typedef std::pair<point, int> value;
	typedef bg::model::box<point> box;

    //[rtree_quickstart_create
    // create the rtree using default constructor
	typedef std::vector<value> container;
	container points;
	
	int n = pts.size();
	
	cout << "point size: " << n << endl;
	
	//vector<point_type> ps = generate_points(n, min_val, max_val);
	//vector<query_type> queries = generate_queries(num_queries, min_val, max_val);
	
	for ( int i = 0 ; i < n ; ++i )
	{
		point p(pts[i].x, pts[i].y);
		points.push_back(std::make_pair(p, i));
	}	
	
    //[rtree_quickstart_insert
    // create some values

	//bgi::rtree< value, bgi::quadratic<16> > rtree(segments.begin(), segments.end());
	
	std::cout << "start building ..." << std::endl;
	parlay::internal::timer t("range", verbose);
	bgi::rtree<value, bgi::linear<16, 4> > rtree(points.begin(), points.end());
	
	cout << "build finished"  << endl;
	t.next("build");
	
	size_t total = 0;
	//int num_queries = queries.size();
	// It takes very long for queries, so using q_num = n/3 may take forever
	int num_queries = 1000;
	size_t* total_array = new size_t[num_queries];
	
	cout << "query size: " << num_queries << endl;
	
	cout << "start querying ..."  << endl;
	#pragma omp parallel for
	for (int i = 0; i < num_queries; i++) {
		//cout << "query: " << queries[i].x1 << " " << queries[i].y1 << " " << queries[i].x2 << " " << queries[i].y2 << endl;
		box query_box(point(queries[i].x1, queries[i].y1), point(queries[i].x2, queries[i].y2));
		std::vector<value> result_s;
		rtree.query(bgi::intersects(query_box), std::back_inserter(result_s));
		//cout << result_s.size() << endl;
		//total += result_s.size();
		total_array[i] = result_s.size();
	}
	//cout << "query time: " << query_tm.stop() << endl;
	t.next("query");
	for (int i = 0; i < num_queries; i++) total+=total_array[i];
	cout << "total: " << total << endl;
    //]
	
    

    // note: in Boost.Geometry WKT representation of a box is polygon

    //[rtree_quickstart_output
    // display results
    //std::cout << "spatial query box:" << std::endl;
    //std::cout << bg::wkt<point>(query_point) << std::endl;
    //std::cout << "spatial query result:" << std::endl;
    //BOOST_FOREACH(value const& v, result_s)
        //std::cout << bg::wkt<box>(v.first) << " - " << v.second << std::endl;

    return total;
}

//]
