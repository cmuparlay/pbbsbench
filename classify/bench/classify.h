using value = unsigned char; 
using values = parlay::sequence<value>;

struct feature {
  bool discrete;
  int num;
  values vals;
  feature(bool discrete, int num) : discrete(discrete), num(num) {}
  feature(bool d, int n, values v) : discrete(d), num(n), vals(v) {}
};

using features = parlay::sequence<feature>;

void classify(features &);
