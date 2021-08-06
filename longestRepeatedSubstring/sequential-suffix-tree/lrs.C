#include <unordered_map>
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

struct SuffixTree {
	const size_t INF = 1e15;	size_t node = 0, pos = 0, cap, sz = 1, n = 0;
	string s;  vector<size_t> len, fpos, link;	vector<unordered_map<size_t, size_t>> to;
	size_t make_node(size_t _pos, size_t _len) {	fpos[sz] = _pos, len [sz] = _len; return sz++; }
	void add_letter(unsigned int c) {
		size_t last = 0;  s += c; n++; pos++;
		while(pos > 0) {
			while(pos > len[to[node][s[n - pos]]]) node=to[node][s[n-pos]], pos-=len[node];
			size_t edge = s[n - pos], &v = to[node][edge], t = s[fpos[v] + pos - 1];
			if (v == 0) v = make_node(n - pos, INF), link[last] = node, last = 0;
			else if (t == c) { link[last] = node;	return;	}
			else {
				size_t u = make_node(fpos[v], pos - 1);
				to[u][c] = make_node(n - 1, INF),	to[u][t] = v;
				fpos[v] += pos - 1,	len [v] -= pos - 1;
				v = u, link[last] = u, last = u;
			}
			if(node == 0) pos--;
			else node = link[node];
		}
	}

	SuffixTree(size_t N) : cap(2*N), len(cap), fpos(cap),
    link(cap), to(cap) { len[0] = INF; s.reserve(N); }
};


// returns
//  1) the length of the longest match
//  2) start of the first string in s
//  3) start of the second string in s
result_type lrs(charseq const &s) {
  parlay::internal::timer t("lrs", true);

  // First, build a suffix tree on the string
  SuffixTree tree(s.size());
  for (const auto& c : s) tree.add_letter(c);
  t.next("build suffix tree");

  // Second, find the deepest internal node of the tree. The path in the tree
  // from the root to this node corresponds to the longest repeated substring
  size_t deepest_internal_node = 0;
  size_t deepest_depth = 0;

  std::function<void(size_t,size_t)> dfs = [&](size_t u, size_t h) {
    bool leaf = true;
    for (const auto& e : tree.to[u]) if (e.second) {
      size_t len = std::min(tree.len[e.second], s.size() - tree.fpos[e.second]);
      dfs(e.second, h + len);
      leaf = false;
    }
    if (!leaf && h > deepest_depth) {
      deepest_internal_node = u;
      deepest_depth = h;
    }
  };
  dfs(0, 0);
  t.next("find deepest internal node");

  // Lastly, check the leaves that are decendants of this deepest internal
  // node, since those correspond to all of the occurences of the repeated
  // substring. Note that every child of the deepest internal node is definitely
  // a leaf, otherwise it wouldn't be the deepest internal node.
  std::vector<size_t> occurences;
  for (const auto& e : tree.to[deepest_internal_node]) if (e.second) {
    size_t position = tree.fpos[e.second] - deepest_depth;
    occurences.push_back(position);
  }
  t.next("find all occurences of repeated substring");
  
  return std::make_tuple(deepest_depth,occurences[0],occurences[1]);
}

