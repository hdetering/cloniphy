#include "amp.hpp"
/*#include <boost/format.hpp>
using boost::format;
using boost::str;*/
#include <cmath>
#include <functional>
//#include <iostream>
using namespace std;

namespace amp {

//! Simulate Multiple Displacement Amplification (MDA).
/*!
  \param cvg output cvg vector of amplicon coverage values (index is position in sequence)
  \param seq_len length of primary sequence
  \param amplicon_size maximum size of generated amplicons
  \param fold factor by which the sequence should be amplified (determines when to stop)
 */
void simulateMda(
  vector<unsigned> &cvg,
  unsigned seq_len,
  double amp_size_mean,
  double amp_size_sd,
  int fold,
  RandomNumberGenerator<> rng
) {
  // coverage evenness (parameter for exponential dist.)
  double lambda = 0.1;
  // number of amplicons needed
  double num_amp = seq_len * fold * lambda / amp_size_mean;
  // per-site probability for amplicon start
  double p = num_amp / double(seq_len);
  function<double()> rf_amp_len = rng.getRandomFunctionGammaMeanSd(amp_size_mean, amp_size_sd);
  function<double()> rf_amp_cvg = rng.getRandomFunctionExponential(lambda);
  function<double()> rf_coin = rng.getRandomIndexWeighted({1-p, p});

  cvg = vector<unsigned>(seq_len, 1); // one copy is the primary sequence itself
  vector<int> v_amps; // active amplicons' remaining lengths
  for (unsigned i=0; i<seq_len; ++i) {
    if (rf_coin()==1) {
      unsigned amp_len = ceil(rf_amp_len());
      int amp_cvg = ceil(rf_amp_cvg());
      for (auto j=0; j<amp_len && j<seq_len; ++j) {
        cvg[i+j] += amp_cvg;
      }
    }
  }
}

void simulateMdaProcessHaploid(
  vector<unsigned> &cvg,
  unsigned seq_len,
  double amp_size_mean,
  double amp_size_sd,
  int fold,
  RandomNumberGenerator<> rng
) {
  int min_amplicon_size = 2000;
  cvg = vector<unsigned>(seq_len, 1); // one copy is the primary sequence itself
  unsigned long total_len = seq_len;
  vector<unsigned> vec_frag_start = { 0 }; // start position of fragments
  vector<unsigned> vec_frag_len = { seq_len }; // length of fragments
  function<double()> rf_amp_len = rng.getRandomFunctionGammaMeanSd(amp_size_mean, amp_size_sd);
  function<short()> rf_strand = rng.getRandomFunctionInt(short(0),short(1));

  while (total_len < seq_len*fold) {
//fprintf(stderr, "%d\n", total_len);
    // pick fragment to amplify
    int idx_frag = rng.getRandomIndexWeighted(vec_frag_len)();
    unsigned src_frag_len = vec_frag_len[idx_frag];
    // pick position in fragment to amplify from
    unsigned pos = rng.getRandomFunctionInt(unsigned(0), src_frag_len-1)();
    // pick strand
    bool is_pstrand = (rf_strand() == 0);

    int new_frag_len = floor(rf_amp_len());
    if (is_pstrand)
      new_frag_len = src_frag_len-pos < new_frag_len ? src_frag_len-pos : new_frag_len;
    else
      new_frag_len = pos+1 < new_frag_len ? pos+1 : new_frag_len;
    if (new_frag_len < min_amplicon_size)
      continue;

    unsigned new_frag_start = vec_frag_start[idx_frag] + pos;
    if (!is_pstrand)
      new_frag_start -= new_frag_len-1;
    // update coverage vector
    for (auto c=cvg.begin()+new_frag_start; c<cvg.begin()+new_frag_start+new_frag_len; ++c) {
      (*c)++;
    }
    // add new amplicon to available fragments
    vec_frag_start.push_back(new_frag_start);
    vec_frag_len.push_back(new_frag_len);
    total_len += new_frag_len;
  }
}

void simulateMdaProcessDiploid(
  vector<vector<unsigned>> &cvg,
  unsigned long long seq_len,
  double amp_size_mean,
  double amp_size_sd,
  int fold,
  RandomNumberGenerator<> rng
) {
  const short ploidy = 2;
  int min_amplicon_size = 2000;
  const unsigned long long target_len = ploidy * seq_len * fold;
  unsigned long long total_len = ploidy*seq_len;
  // two copies of primary sequence are present at start
  cvg = vector<vector<unsigned>>(ploidy, vector<unsigned>(seq_len, 1));
  // start position of fragments
  vector<unsigned> vec_frag_start = vector<unsigned>(ploidy, 0);
  // fragment lengths
  vector<unsigned> vec_frag_len = vector<unsigned>(ploidy, seq_len);
  vector<unsigned short> vec_frag_allele = { 0, 1 }; // allelic source (\in {0,1})
  function<double()> rf_amp_len = rng.getRandomFunctionGammaMeanSd(amp_size_mean, amp_size_sd);
  function<short()> rf_strand = rng.getRandomFunctionInt(short(0),short(1));

  while (total_len < target_len) {
//fprintf(stderr, "%d\n", total_len);
    // pick fragment to amplify
    int idx_frag = rng.getRandomIndexWeighted(vec_frag_len)();
    unsigned src_frag_len = vec_frag_len[idx_frag];
    unsigned short src_frag_allele = vec_frag_allele[idx_frag];
    // pick position in fragment to amplify from
    unsigned pos = rng.getRandomFunctionInt(unsigned(0), src_frag_len-1)();
    // pick strand
    bool is_pstrand = (rf_strand() == 0);

    int new_frag_len = floor(rf_amp_len());
    if (is_pstrand)
      new_frag_len = src_frag_len-pos < new_frag_len ? src_frag_len-pos : new_frag_len;
    else
      new_frag_len = pos+1 < new_frag_len ? pos+1 : new_frag_len;
    if (new_frag_len < min_amplicon_size)
      continue;

    unsigned new_frag_start = vec_frag_start[idx_frag] + pos;
    if (!is_pstrand)
      new_frag_start -= new_frag_len-1;
    // update coverage vector
    for (auto it_cvg = cvg[src_frag_allele].begin()+new_frag_start;
         it_cvg < cvg[src_frag_allele].begin()+new_frag_start+new_frag_len;
         ++it_cvg) {
      (*it_cvg)++;
    }
    // add new amplicon to available fragments
    vec_frag_start.push_back(new_frag_start);
    vec_frag_len.push_back(new_frag_len);
    vec_frag_allele.push_back(src_frag_allele);
    total_len += new_frag_len;
  }
}


} /* namespace amp */
