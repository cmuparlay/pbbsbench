/*****************************************************************************/
/*                                                                           */
/*  Header: fst.hpp                                                          */
/*                                                                           */
/*  Accompanies STANN Version 0.74                                           */
/*  Nov 15, 2010                                                             */
/*                                                                           */
/*  Copyright 2007, 2008                                                     */
/*  Michael Connor and Piyush Kumar                                          */
/*  Florida State University                                                 */
/*  Tallahassee FL, 32306-4532                                               */
/*                                                                           */
/*****************************************************************************/

#ifndef __STANN_FAIR_SPLT_TREE__
#define __STANN_FAIR_SPLT_TREE__

#include<limits>
#include<vector>
#include<iostream>

template<typename Point>
class dim_sort_pred
{
public:
  int dim;
  dim_sort_pred(int D)
  {
    dim=D;
  }
  bool operator()(const Point &a, const Point &b)
  {
    return a[dim]< b[dim];
  }
};

template<typename Point>
class fst_node
{
  typedef typename std::vector<fst_node<Point> >::size_type size_type;
  typedef typename std::vector<Point>::iterator PItr;

public:
  
  Point lc;  //Lower corner of the bounding box
  Point uc;  //Upper corner of the bounding box
  PItr first; //Smallest point, according to the splitting dimenstion
  PItr last; //1 past the last point, according to splitting dimension
  PItr cut; //Point at which the node is split
  double radius; //Radius of bounding box.  Equal to 1/2 the squared corner to corner distance
  
  fst_node<Point>* left; //left child of the node in the FST
  fst_node<Point>* right; //right child of the node in the FST
  int size;

  double node_distance(const fst_node<Point> *b)
  {
    double dist=0;
    for(unsigned int j=0;j < Point::__DIM;++j)
      {
	if(uc[j] < b->lc[j])
	  dist += (b->lc[j]-uc[j])*(b->lc[j]-uc[j]);
	else if(b->uc[j] < lc[j])
	  dist += (lc[j]-b->uc[j])*(lc[j]-b->uc[j]);
      }
    return dist;
  }
};
template<typename Point>
class fair_split_tree
{

  typedef typename std::vector<fst_node<Point> >::size_type size_type;
  typedef typename std::vector<Point>::iterator PItr;
  
public:
  
  fst_node<Point> root;
  int DIM; 

  fair_split_tree(std::vector<Point> &points)
  {
    DIM = Point::__DIM;
    build_tree(points.begin(), points.end(), &root);
  }
  
  void build_tree(PItr begin, PItr end, fst_node<Point> *node)
  {
    node->first=begin;
    node->last=end;
    node->size = (int) (end-begin);
    if(node->size==0)
      {
	std::cout << "Error, no size!" << std::endl;
	return;
      }
    //if we're in an internal node
    if(node->size > 1)
      {
	int spread_dim;

	//Find the max and min coordinates
	PItr I = begin;
	for(int i=0;i < DIM;++i)
	  {
	    node->lc[i] = (*I)[i];
	    node->uc[i] = (*I)[i];
	  }
	I++;
	for(;I != end;++I)
	  {
	    for(int i=0;i < DIM;++i)
	      {
		if(node->lc[i] > (*I)[i])
		  node->lc[i] = (*I)[i];
		if(node->uc[i] < (*I)[i])
		  node->uc[i] = (*I)[i];
	      }
	  }
	//Calculate the radius of the node
	node->radius=0;
	for(int i=0;i < DIM;++i)
	{
	  node->radius= node->lc.sqr_dist(node->uc);
	}
	node->radius = node->radius*0.5;
	
	//std::cout << "LC: " << node->lc << std::endl;
	//std::cout << "UC: " << node->uc << std::endl;
	//std::cout << node->radius << std::endl;
	//exit(0);
	//Determine the dimension of greatest spread
	double spread;
	spread = node->uc[0]-node->lc[0];
	spread_dim=0;
	for(int i=1;i < DIM;++i)
	  {
	    if((node->uc[i]-node->lc[i])> spread)
	      {
		spread = node->uc[i]-node->lc[i];
		spread_dim=i;
	      }
	  }

	//sort the points in the dimension of greatest spread
	dim_sort_pred<Point> lt(spread_dim);
	sort(begin, end, lt);

	//determine the cut
	spread=(spread*0.5)+(node->lc[spread_dim]);
	if(node->size == 2)
	  {
	    node->cut = begin+1;
	  }
	else
	  {
	    node->cut = begin+(node->size/2);
	    if((*node->cut)[spread_dim] < spread)
	      while((*node->cut)[spread_dim] < spread)
		++node->cut;
	    else
	      {
		while((*node->cut)[spread_dim] > spread)
		  --node->cut;
		++node->cut;
	      }
	  }
	//if(node->cut == end) std::cout << "Error, cut=end" << std::endl;
	//if(node->cut == begin) std::cout << "Error, cut=begin" << std::endl;
	//Recurse
	node->left = new fst_node<Point>;
	node->right = new fst_node<Point>;
	build_tree(begin, node->cut, node->left);
	build_tree(node->cut, end, node->right);
      }
    //Otherwise, we're in a leaf
    else
      {
	node->radius=0;
	node->uc = (*begin);
	node->lc = (*begin);
	node->left = NULL;
	node->right=NULL;
      }
  }
  
};
#endif
