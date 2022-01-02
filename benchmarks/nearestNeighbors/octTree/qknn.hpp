#include <vector>
#include <queue>
#include <algorithm>

/*! \file qknn.hpp
\brief Implements priority queue functions for dpoints and point pairs */

//! Priority Queue element comparator
/*! This class orders priority queue elements based on the distance
  given as the first item in the pair 
*/

template <class vtx>
class q_intelementCompare {   
public:

  //! Less than operator
  /*! 
    Compares two priority queue elements based on thier distance
    \param p1 First element to be compared
    \param p2 Second element to be compared
    \return Returns true if p1 distance is less than p2 distance
  */
  bool operator()( const std::pair<vtx*, double> p1,  
                   const std::pair<vtx*, double> p2 ){
    return p1.second < p2.second;
  }
};

//! Distance Priority Queue
/*! 
  Implements a priority queue for pairs of floating point 
  distances and array indexes.  The priority queue is ordered 
  based on the squared distance stored in the first element
  of the pair.
*/

template <class vtx>
class qknn 
{
private:
  long unsigned int K;
  typedef std::pair<vtx*, double> q_intelement;
  //typedef std::priority_queue<q_intelement, parlay::sequence<q_intelement>, q_intelementCompare>
  typedef std::priority_queue<q_intelement, std::vector<q_intelement>, q_intelementCompare<vtx>>  
  PQ;
  PQ pq;
  
public:
  
  //! Constructor
  /*! 
    Creates an empty priority  queue.
   */
  qknn(){};
  
  //! Largest distance
  /*! 
    Returns the largest distance value stored in the priority queue
    \return Largest distance value
  */

  void pop(){
    pq.pop();
  }

  q_intelement top(){
    return pq.top();
  }

  double topdist(void)
  {
    return pq.top().second;
  }
  
  //! Set Size
  /*! 
    Sets the size of the priority queue.  This should be set before the
    queue is used
    \param k The maximum number of elements to be stored in the queue.
  */
  void set_size(long unsigned int k)
  {
    K = k;
  }
  
  /*bool is_element(double dist, long int p)
  {
  }*/
  //! Point with largest distance
  /*!
    Returns the index associated with the largest element in the queue.
    \return Index of largest (most distant) element
  */



  //! Update queue
  /*! 
    Updates the queue with the given distance and point
    \param dist Distance of point to be added
    \param p index of point to be added
    \return True if a point was added to the queue
  */
  bool update(vtx* v, double d)
  {
    if(size() < K)
      {
	q_intelement tq(v, d);
	pq.push(tq);
	return true;
      }
    else if(topdist() > d)
      {
	pq.pop();
	q_intelement tq(v, d);
	pq.push(tq);
	return true;
      }
    return false;
  }


  //! Size function
  /*! 
    Returns the current size of the queue 
    \return Size
  */  
   long unsigned int size(){ return pq.size(); }
};

