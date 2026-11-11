// Paired TCRdist kernel, a line-by-line port of pwseqdist's
// `_nb_vector_tcrdist` (used by tcrdist3, MIT License,
// Copyright (c) 2020 Koshlan Mayer-Blackwell, Andrew Fiore-Gartland).
//
// Sequences are encoded once (encode_seqs) as 0..19 for the 20 amino acids
// and 20 for anything else (gap '.', 'X', ...). As in tcrdist3, the distance
// between "anything else" and any symbol is 0.
#include <Rcpp.h>
#include <string>
#include <vector>
#include <algorithm>
using namespace Rcpp;

static const char *AA = "ARNDCQEGHILKMFPSTWYV";

static int aa_index(char c) {
  for (int i = 0; i < 20; i++) if (AA[i] == c) return i;
  return 20;
}

// [[Rcpp::export(.encode_seqs)]]
List encode_seqs(CharacterVector seqs) {
  int n = seqs.size();
  IntegerVector offs(n + 1);
  IntegerVector lens(n);
  std::vector<int> codes;
  codes.reserve(n * 16);
  offs[0] = 0;
  for (int i = 0; i < n; i++) {
    if (CharacterVector::is_na(seqs[i])) {
      lens[i] = NA_INTEGER;
      offs[i + 1] = offs[i];
      continue;
    }
    std::string s = Rcpp::as<std::string>(seqs[i]);
    lens[i] = s.size();
    for (size_t j = 0; j < s.size(); j++) codes.push_back(aa_index(s[j]));
    offs[i + 1] = offs[i] + s.size();
  }
  return List::create(_["codes"] = IntegerVector(codes.begin(), codes.end()),
                      _["offsets"] = offs, _["lengths"] = lens);
}

// Distance from one query sequence to every encoded sequence.
// dm: 21 x 21 integer distance matrix (column-major, last row/col zero).
// Returns min_dist * dist_weight + len_diff * gap_penalty (pwseqdist), NA for NA seqs.
// [[Rcpp::export(.tcrdist_one_to_many)]]
IntegerVector tcrdist_one_to_many(IntegerVector query, IntegerVector codes,
                                  IntegerVector offsets, IntegerVector lengths,
                                  IntegerMatrix dm, int dist_weight,
                                  int gap_penalty, int ntrim, int ctrim,
                                  bool fixed_gappos) {
  int n = lengths.size();
  IntegerVector out(n);
  int q_L = query.size();
  const int *q = query.begin();
  const int *cd = codes.begin();
  const int *D = dm.begin();
  const int nr = dm.nrow();
  for (int k = 0; k < n; k++) {
    int s_L = lengths[k];
    if (s_L == NA_INTEGER) { out[k] = NA_INTEGER; continue; }
    const int *s = cd + offsets[k];
    if (q_L == s_L) {
      int d = 0;
      for (int i = ntrim; i < q_L - ctrim; i++) d += D[q[i] + nr * s[i]] * dist_weight;
      out[k] = d;
      continue;
    }
    int short_len = std::min(q_L, s_L);
    int len_diff = std::abs(q_L - s_L);
    int min_gappos, max_gappos;
    if (fixed_gappos) {
      // Python floor division: (short_len - 5) // 2
      int t = short_len - 5;
      int fl = (t >= 0) ? t / 2 : -((-t + 1) / 2);
      min_gappos = std::min(6, 3 + fl);
      max_gappos = min_gappos;
    } else {
      min_gappos = 5;
      max_gappos = short_len - 1 - 4;
      while (min_gappos > max_gappos) { min_gappos -= 1; max_gappos += 1; }
    }
    int min_dist = -1;
    for (int gappos = min_gappos; gappos <= max_gappos; gappos++) {
      int tmp = 0;
      int remainder = short_len - gappos;
      for (int n_i = ntrim; n_i < gappos; n_i++) tmp += D[q[n_i] + nr * s[n_i]];
      for (int c_i = ctrim; c_i < remainder; c_i++)
        tmp += D[q[q_L - 1 - c_i] + nr * s[s_L - 1 - c_i]];
      if (tmp < min_dist || min_dist == -1) min_dist = tmp;
      if (min_dist == 0) break;
    }
    out[k] = min_dist * dist_weight + len_diff * gap_penalty;
  }
  return out;
}
