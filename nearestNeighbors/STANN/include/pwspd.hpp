/*****************************************************************************/
/*                                                                           */
/*  Header: pwspd.hpp                                                        */
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

#ifndef __STANN_PWSPD__
#define __STANN_PWSPD__

#include <cstring>
#include <limits>
#include <vector>
#include <fst.hpp>
#include <zorder_lt.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

/*! \file
  \brief Implementation of a parallel well seperated pair decomposition for a compressed quad tree
*/

/*! Well Seperated Pair
  \brief Implements functions needed by a Well Seperated Pair object
*/
template<typename Point>
class WSP
{
  typedef typename std::vector<Point>::size_type size_type;
  typedef typename std::vector<Point>::iterator PItr;
  typedef fst_node<Point>* index;

public:

  //Before the BCCP is computed, first and last point
  //To the quad-tree boxes that make up the WSP
  //After BCCP is computed, they point to the two points
  //that make up the BCCP
  /*! first
    \brief Indicates either the first quadtree box, or point
    Initially, this pointer indcates the comp_quad_tree_node that
    is the first half of the well-seperated pair.  After the BCCP is
    calculated, this because a pointer to a point, indicating it is
    half of the BCCP
  */
  void* first;
  /*! second
    \brief Indicates either the second quadtree box, or point
    Initially, this pointer indcates the comp_quad_tree_node that
    is the second half of the well-seperated pair.  After the BCCP is
    calculated, this because a pointer to a point, indicating it is
    half of the BCCP
  */
  void* second;
  /*! bccp_dist
    \brief Indicates either the BCCP distance, or the marginal BCCP distance
    If this value is negative, it indicates the current minimum possible BCCP distance known.  If it is positive, it indicates the actual BCCP distance
  */
  double bccp_dist;
  /*! Constructor
    \brief Default constructor
  */
  WSP()
  {    bccp_dist = -1;
  };

  /*! Constructor
    \brief Construcs a WSP consisting of two comp_quad_tree_node elements
    \param f First comp_quad_tree_node
    \param s Second comp_quad_tree_node
  */
  WSP(index f, index s) : first(f), second(s)
  {

    //We store the marginal distance as a negative to indicate
    //That we haven't computed the actual BCCP yet
    bccp_dist = -f->node_distance(s);
    if(bccp_dist == 0)
      {
	std::cout << std::endl << "Error in WSP construction!" << std::endl;
	//std::cout << (*(f->lc))[0] << "," << (*(f->low))[1]  << std::endl;
	//std::cout << (*(s->lc))[0] << "," << (*(s->low))[1] << std::endl;
	if(f==s)
	  std::cout << "F==S!" << std::endl;
	exit(0);
      }
  };
  
  /*! size
    \brief returns the number of points in the WSP
  */
  long int size()
  {
    if(bccp_done())
      return 0;
    //std::cout << "bccp_dist: " << bccp_dist << std::endl;
    return ((index)first)->size+((index)second)->size;
  };

  /*! get bccp distance
    \brief returns the currently computed BCCP distance
  */
  double get_bccp_dist() const
  {
    return (double) fabs(bccp_dist);
  }

  /*! BCCP done
    \brief returns true if the actual BCCP has been computed, false otherwise
  */
  bool bccp_done() const
  {
    return bccp_dist >= 0;
  }

  /*! Assignment Operator
    \brief Assignment operator
  */
  WSP<Point>& operator=(const WSP<Point> b)
  {
    first = b.first;
    second = b.second;
    bccp_dist = b.bccp_dist;
    return *this;
  }

  /*! BCCP
    \brief computes the BCCP if it has not been
    For WSPs of small size (<32) this computes the BCCP by brute force.
    Otherwise it calls a recursive algorithm.
  */
  void bccp()
  {
    if(bccp_done())
      return;    
    //Note: After BCCP is computed, the void* first and second
    //point to Points, not comp_quad_tree_nodes
    
    //For small sized quad-tree boxes, 
    //Brute force computation of BCCP is faster
    if(size() < 32)
      {
	brute_force_bccp((index)first, (index)second);
	return;
      }

    //For small sized boxes, compute using 
    //recursive algorithm
    do_bccp();
    return;
  };
private:
  //Compute the BCCP by brute force
  //Note: This works well because
  //the points in the quad-tree boxes
  //are already in a contiguous chunk the the
  //z-order sorted point list
  void brute_force_bccp(index l, index r)
  {
    double dist;
    bccp_dist = (std::numeric_limits<double>::max)();
    for(PItr i=l->first;i != l->last;i+=1)
      {
	for(PItr j=r->first;j != r->last;j+=1)
	  {
	    dist = i->sqr_dist(*j);
	    if(dist < bccp_dist)
	      {
		bccp_dist=dist;
		first=&*i;
		second=&*j;
	      }
	  }
      }
  };
  
  //Recursive BCCP computation.
  //Decends the compressed quadtree and 
  //computes the bccp
  void bccp_recurse(index r, index l)
  {
    double d1, d2;
    index t;
    if((r->radius==0) && (l->radius==0))
      {
	d1 = r->node_distance(l);
	if(d1 < bccp_dist)
	  {
	    bccp_dist=d1;
	    first=&*(r->first);
	    second=&*(l->first);
	  }
	return;
      }
    if(r->radius < l->radius)
      {
	t=l;
	l=r;
	r=t;
      }

    d1 = l->node_distance(r->left);
    d2 = l->node_distance(r->right);

    if(d1 < d2)
      {
	if(d1 < bccp_dist)
	  {
	    bccp_recurse(r->left, l);
	  }
	if(d2 < bccp_dist)
	  {
	    bccp_recurse(r->right, l);
	  }
      }
    else
      {
	if(d2 < bccp_dist)
	  {
	    bccp_recurse(r->right, l);
	  }
	if(d1 < bccp_dist)
	  {
	    bccp_recurse(r->left, l);
	  }
      }
  }

  //compute the bi-chromatic closest pair from two well seperated ranges
  void do_bccp()
  {
    bccp_dist = ((index)first)->first->sqr_dist(*(((index)second)->first));
    index f,s;
    f = (index) first;
    s = (index) second;
    first = &*(f->first);
    second = &*(s->first);
    bccp_recurse(f,s);
  }
};

/*! Parallel Well Seperated Pair Decomposition Class
  \brief Implements a prallel well-seperated pair decomposition of a point set
*/
template<typename Point>
class pwspd
{
public:
  typedef typename std::vector<Point>::size_type size_type;
  typedef fst_node<Point>* Node;
  typedef typename std::vector<Point>::iterator PItr;
  typedef std::pair<PItr, PItr> Edge;

  std::vector<Point> &Points;
  fair_split_tree<Point> &cqtree;
  double spread;
  /*! Constructor
    \brief Constructor
    \param P vector of points for which the WSPD is to be constructed.
    \param S seperation factor for the WSPD (boxes will be considered well separated if they have a minimum distance greater than S times the maximum corner to corner distance of the boxes)
  */
  pwspd(std::vector<Point> &P, fair_split_tree<Point> &T, double S=2.0) : Points(P), cqtree(T)
  {
    setSpread(S);
  }
  /*! Destructor
    \brief Destructor
  */
  ~pwspd()
  {
  };
  
  /*! Set Spread
    \brief sets the spread of the WSPD
    \param S the spread factor
    Note, this will only effect the WSPD if set before the run command is issued
  */
  void setSpread(double S)
  {
    spread=S;
  };

  /* Get Spread
     \brief returns the spread of the WSPD
  */
  double getSpread()
  {
    return sqrt(spread);
  };
  
  /*! Run
    \brief Generates the WSPD
    \param Wsp Vector of WSPs storing the output
    \param num_threads Number of threads to use (requires OPENMP)
  */
  void run(std::vector<WSP<Point> > &Wsp, int num_threads=1)
  {
    recurse(Wsp, &(cqtree.root), num_threads);
  };

  void recurse(std::vector<WSP<Point> > &Wsp, fst_node<Point> *node, int num_threads)
  {
    if(node->left->size > 1)
      recurse(Wsp, node->left, num_threads);
    if(node->right->size > 1)
      recurse(Wsp, node->right, num_threads);
    decompose(node->left, node->right, Wsp);
}
  void decompose(fst_node<Point> *left, fst_node<Point> *right, std::vector<WSP<Point> > &Wsp)
  {
    double dist = left->node_distance(right);
    double radius = left->radius;

    if(right->radius > radius)
      radius = right->radius;
    
    if(dist > (spread*radius))
      {
	Wsp.push_back(WSP<Point>(left, right));
      }
    else if(left->radius > right->radius)
      {
	decompose(left->left, right, Wsp);
	decompose(left->right, right, Wsp);
      }
    else
      {
	decompose(right->left, left, Wsp);
	decompose(right->right, left, Wsp);
      }
      
  }
};
#endif
