#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "parlay/sequence.h"
#include "parlay/internal/get_time.h"

using std::string;
using std::vector;
using std::unordered_map;

using charseq = parlay::sequence<unsigned char>;
using result_type = std::tuple<size_t,size_t,size_t>;

// Suffix array construction with LCP. The suffix array is built into sarray and the
// LCP into lcp. NOTE: sarray does not include the empty suffix. lcp[i] is the longest
// common prefix between the strings at sarray[i-1] and sarray[i], lcp[0] = 0.
// Complexity: O(N) or O(N log(N)) for suffix array. O(N) for LCP.

template<typename STRING, typename int_t>
struct suffix_array {
  const STRING& str;
  int_t n; vector<int_t> sarray, lcp;
  
  void bucket(vector<int_t>& a, vector<int_t>& b, vector<int_t>& r, int_t n, int_t K, int_t off=0) {
    vector<int_t> c(K+1, 0);
    for (int_t i=0; i<n; i++) c[r[a[i]+off]]++;
    for (int_t i=0, sum=0; i<=K; i++) { int_t t = c[i]; c[i] = sum; sum += t; }
    for (int_t i=0; i<n; i++) b[c[r[a[i]+off]]++] = a[i];
  }

  suffix_array(const STRING& s) : n(s.size()), str(s) { build_sarray(); build_lcp(); }
  
  using tiii = std::tuple<int_t,int_t,int_t>;
  using pii = std::pair<int_t,int_t>;
  
  void sarray_int(vector<int_t> &s, vector<int_t> &SA, int_t n, int_t K) {
    int_t n0=(n+2)/3, n1=(n+1)/3, n2=n/3, n02=n0+n2, name=0, c0=-1, c1=-1, c2=-1;
    vector<int_t> s12(n02 + 3, 0), SA12(n02 + 3, 0), s0(n0), SA0(n0);
    for (int_t i=0, j=0; i < n+(n0-n1); i++) if (i%3 != 0) s12[j++] = i;
    bucket(s12, SA12, s, n02, K, 2), bucket(SA12, s12, s, n02, K, 1);
    bucket(s12, SA12, s, n02, K, 0);
    for (int_t i = 0; i < n02; i++) {
      if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2) 
        name++, c0 = s[SA12[i]], c1 = s[SA12[i]+1], c2 = s[SA12[i]+2];
      if (SA12[i] % 3 == 1) s12[SA12[i]/3] = name;
      else s12[SA12[i]/3 + n0] = name;
    }
    if (name < n02) {
      sarray_int(s12, SA12, n02, name);
      for (int_t i = 0; i < n02; i++) s12[SA12[i]] = i + 1;
    } else for (int_t i = 0; i < n02; i++) SA12[s12[i] - 1] = i;
    for (int_t i=0, j=0; i < n02; i++) if (SA12[i] < n0) s0[j++] = 3*SA12[i];
    bucket(s0, SA0, s, n0, K);
    for (int_t p=0, t=n0-n1, k=0; k < n; k++) {
      int_t i = (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2), j = SA0[p];
      if (SA12[t] < n0 ?
          (pii(s[i], s12[SA12[t] + n0]) < pii(s[j], s12[j/3])) :
          (tiii(s[i],s[i+1],s12[SA12[t]-n0+1]) < tiii(s[j],s[j+1],s12[j/3+n0]))) {
        SA[k] = i; t++;
        if (t == n02) for (k++; p < n0; p++, k++) SA[k] = SA0[p];
      } else {
        SA[k] = j; p++;
        if (p == n0) for (k++; t < n02; t++, k++) SA[k] = (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2);
      }
    }
  }
  void build_sarray() {
    if (n <= 1) { sarray.assign(n, 0); return; }
    vector<int_t> s(n+3, 0); sarray.assign(n+3, 0);
    for (int_t i=0; i<n; i++) s[i] = (int_t)str[i] - CHAR_MIN + 1;
    sarray_int(s, sarray, n, 256), sarray.resize(n);
  }

  void build_lcp() {  
    int_t h = 0; vector<int_t> rank(n); lcp.assign(n, 0);
    for (int_t i = 0; i < n; i++) rank[sarray[i]] = i;
    for (int_t i = 0; i < n; i++) {
      if (rank[i] > 0) {
        int_t j = sarray[rank[i]-1];
        while (i + h < n && j + h < n && str[i+h] == str[j+h]) h++;
        lcp[rank[i]] = h;
      }
      if (h > 0) h--;
    }
  } 
};


// returns
//  1) the length of the longest match
//  2) start of the first string in s
//  3) start of the second string in s
template <typename int_t>
result_type lrs_(charseq const &s) {
  parlay::internal::timer t("lrs", true);

  std::cout << "n = " << s.size() << std::endl;

  // First, build a suffix tree on the string
  suffix_array<charseq,int_t> sa(s);
  sa.build_sarray();
  t.next("build suffix array");

  sa.build_lcp();
  t.next("build LCP");
  
  size_t idx = std::distance(sa.lcp.begin(), std::max_element(sa.lcp.begin(), sa.lcp.end()));
  t.next("max element");
    
  return std::make_tuple(sa.lcp[idx], sa.sarray[idx], sa.sarray[idx+1]);
}

result_type lrs(charseq const &s) {
  return lrs_<unsigned int>(s);
}

