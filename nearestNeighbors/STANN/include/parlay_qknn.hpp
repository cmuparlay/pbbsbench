class qknn {
private:
  size_t k;
  int num_elts;
  typedef std::pair<double, long int> q_intelement;
  parlay::sequence<q_intelement> pq;
  double epsilon;
  
public:
  
  qknn(size_t k, double epsilon=1.0) :
    k(k), num_elts(0), pq(parlay::sequence<q_intelement>(k)), epsilon(epsilon) {};
  
  double topdist(void) {
    return pq[0].first;
  }
  
  long int top() {
    return pq[0].second;
  }

  bool update(double dist, long int p) {
    if (num_elts < k) {
      pq[num_elts++] = q_intelement(dist, p);
      int i = num_elts-1;
      while (i > 0 && pq[i-1].first < pq[i].first) {
	swap(pq[i-1], pq[i]);
	i--;
      }
      return true;
    } else if (topdist()*epsilon > dist) {
      pq[0] = q_intelement(dist, p);
      int i = 1;
      while (i < num_elts && pq[i-1].first < pq[i].first) {
	swap(pq[i-1], pq[i]);
	i++;
      }
      return true;
    } else return false;
  }
  
  template <typename Vect>
  void answer(Vect& pl) {
    pl.resize(k);
    for (int i=0; i < num_elts; i++)
      pl[i] = pq[num_elts-i-1].second;
  };

  void answer(std::vector<long unsigned int>& pl, std::vector<double> &pd) {    // not supported
    abort();
  }

  long unsigned int size() { return num_elts; }
};
