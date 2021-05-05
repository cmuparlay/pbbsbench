#include <iostream>
#include <algorithm>
#include <cctype>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/io.h"
#include "parlay/random.h"
#include "common/parse_command_line.h"

int max_continuous = 64;
int max_discrete = 256;
using parlay::map;
using parlay::reduce;
using parlay::tabulate;
using parlay::sequence;

template <typename Seq, typename F>
bool all(Seq const &s, F const &f) {
  return (parlay::count_if(s, f) == s.size());
}

using charseq = sequence<char>;

template <typename Seq>
sequence<sequence<charseq>> read_csv(Seq const &str) {
  auto is_line = [] (char c) -> bool { return c == '\r' || c == '\n'|| c == 0;};
  auto is_item = [] (char c) -> bool { return c == ',';};
  auto is_num = [] (sequence<char> str) -> bool { 
    return all(str, [] (char c) {return std::isdigit(c) || c=='-' || c=='.';});};
  auto is_int = [] (sequence<char> str) -> bool { 
    return all(str, [] (char c) {return std::isdigit(c) || c=='-';});};

  auto process_line = [&] (auto line) {
      return parlay::tokens(line, is_item);};

  auto x = parlay::map_tokens(str, process_line, is_line);
  int num_columns = x[0].size();

  if (!all(x, [&] (auto const &s) {return s.size() == num_columns;})) {
    std::cout << "rows are not all the same size" << endl;
    return sequence<sequence<charseq>>();
  }
  return x;
}

auto convert(sequence<sequence<charseq>> const &rows) {
  auto header = rows[0];
  size_t num_features = header.size();
  size_t num_entries = rows.size() - 1;

  // transpose
  auto features = tabulate(num_features, [&] (size_t i) {
                    return tabulate(num_entries, [&] (size_t j) {
                      return std::move(rows[j+1][i]);});});
   
  auto process_feature = [&] (size_t i) {
      char type = header[i][0];
      auto &feature_vector = features[i];
      auto adjust_range = [&] (auto vals, bool is_discrete=false) {
	auto sorted =  parlay::sort(vals);
	auto uni = parlay::unique(sorted);
	size_t num = uni.size();
	if (is_discrete && num >= max_discrete) {
	  std::cout << "too many discrete values for a feature at index " << i << ", "
		    << num << " found" << std::endl;
	  return sequence<int>(num_entries, 0);
	} else {
	  return map(vals, [&] (auto v) -> int {
	      int j = std::lower_bound(uni.begin(), uni.end(), v) - uni.begin();
	      return (num < max_continuous) ? j : (((long) j * max_continuous) / num);});
	}
      };
      if (type == 'i') { // integers
	return adjust_range(map(feature_vector, parlay::chars_to_int));
      } else if (type == 'f') { // floats
	return adjust_range(map(feature_vector, parlay::chars_to_float));
      } else if (type == 'd') { // discrete char strings
	return adjust_range(feature_vector, true);
      } else {
	std::cout << "invalid type at index " << i << std::endl;
	return sequence<int>(num_entries, 0);
      }
  };

  auto z = tabulate(num_features, process_feature);

  // transpose back
  auto entries = tabulate(num_entries, [&] (size_t i) {
                    return tabulate(num_features, [&] (size_t j) {
                      return std::move(z[j][i]);});});

  return std::pair(header, entries);
}

auto strip_target_label(sequence<sequence<int>> const &rows) {
  return std::pair(map(rows, [&] (sequence<int> const &r) {return to_sequence(r.cut(0,r.size()-1));}),
		   map(rows, [&] (sequence<int> const &r) {return r[r.size()-1];}));
}

sequence<char> comma(1, ',');
sequence<char> newline(1, '\n');

template <typename Seq>
sequence<char> csv_row(Seq const &e) {
  size_t len = e.size()*2;
  auto ss = tabulate(len, [&] (size_t i) {
    if (i == len-1) return newline;
    if (i%2 == 1) return comma;
    return parlay::to_chars(e[i/2]);});
  return parlay::flatten(ss);
};

template <typename Seq>
sequence<char> csv_string(sequence<charseq> const &header,
			    Seq const &rows) {
  return parlay::append(csv_row(header), 
			parlay::flatten(map(rows, [] (auto x) {return csv_row(x);})));
}

int main(int argc, char* argv[]) {
  commandLine P(argc,argv, "[-f <testpercent>] <inFile> <outFile>");
  std::pair<char*,char*> fnames = P.IOFileNames();
  std::string iFile = fnames.first;
  std::string oFile = fnames.second;
  int test_percent = P.getOptionIntValue("-f",20);

  auto str = parlay::file_map(iFile);
  //if (str.size() == 0) {
  //  std::cout << "no such input file, or is empty" << endl;
  //  return 1;
  //}

  auto [header, rows] = convert(read_csv(str));
  size_t n = rows.size();
  auto perm = parlay::random_permutation(n);
  size_t num_train = round(n * (1.0 - test_percent/100.0));
  auto train = map(perm.cut(0,num_train), [&] (size_t i) {return rows[i];});
  auto testr = map(perm.cut(num_train, n), [&] (size_t i) {return rows[i];});
  auto [test, labels] = strip_target_label(testr);

  std::string train_file = oFile;
  std::string test_file = oFile;
  std::string label_file = oFile;
  chars_to_file(csv_string(header, train), train_file.append(".data.train"));
  chars_to_file(csv_string(header, test), test_file.append(".data.test"));
  chars_to_file(csv_row(labels), label_file.append(".data.labels"));
}
