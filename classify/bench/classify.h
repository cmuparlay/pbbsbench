constexpr int max_value = 255;
using value = unsigned char; 
using row = parlay::sequence<value>;
using rows = parlay::sequence<row>;

struct feature {
  bool discrete; // discrete (true) or continuous (false)
  int num;       // max value of feature
  row vals;      // the sequence of values for the feature
  feature(bool discrete, int num) : discrete(discrete), num(num) {}
  feature(bool d, int n, row v) : discrete(d), num(n), vals(v) {}
};

using features = parlay::sequence<feature>;

row classify(features const &Train, rows const &Test, bool verbose);
