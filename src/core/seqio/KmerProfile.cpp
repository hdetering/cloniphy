#include "KmerProfile.hpp"
#include <cmath>

using namespace std;

namespace seqio {

KmerProfile::KmerProfile() {}

KmerProfile::KmerProfile ( const string fn_csv ) {
  this->readCSV(fn_csv);
  this->indexPrefixes();
}

bool 
KmerProfile::readCSV (
  const string filename
) {
  double sum_freqs = 0.0;
  
  // READ DATA
  //---------------------------------------------------------------------------
  
  // read frequency profile from file
  vector<vector<string>> csv_data;
  stringio::readCSV(csv_data, filename, ',');

  bool is_first_kmer = true;
  for (auto row : csv_data) {
    string kmer = row[0];
    double freq;
    // get frequency from string
    try {
      freq = stod(row[1]);
    } catch (invalid_argument& e) {
      fprintf(stderr, "[ERROR] Failed to parse kmer profile entry:\n");
      fprintf(stderr, "\tkey:   '%s'\n", row[0].c_str());
      fprintf(stderr, "\tvalue: '%s'\n", row[1].c_str());
      return false;
    }
    // store kmer length
    if ( is_first_kmer ) {
      kmer_length = kmer.length();
      is_first_kmer = false;
    }

    num_kmers++;
    sum_freqs += freq;
    m_vec_kmers.push_back(kmer);
    m_vec_weight.push_back(freq);
  }

  // SANITY CHECKS
  //---------------------------------------------------------------------------
  // do kmer frequencies contain data?
  if ( num_kmers == 0) {
    fprintf(stderr, "[ERROR] kmer frequencies are empty.\n");
    return false;
  }
  // do all kmers have the same length?
  for (unsigned i=0; i<num_kmers; i++) {
    if ( m_vec_kmers[i].length() != kmer_length ) {
      fprintf(stderr, "[ERROR] kmers are expected to have same length.\n");
      return false;
    }
  }
  // does number of kmers indicate that all kmers are present?
  unsigned num_kmers_exp = pow(4, kmer_length);
  if ( num_kmers != num_kmers_exp ) {
    fprintf(stderr, "[ERROR] number of kmers (%u) not as expected (%d).\n", num_kmers, num_kmers_exp);
    return false;
  }

  // normalize frequency values to sum up to 1.
  // TODO: is this really necessary?
  // fprintf(stderr, "[INFO] All good with kmer spectrum. Converting frequencies to probabilities.\n");
  // for (unsigned i=0; i<num_kmers; i++) {
  //   m_vec_weight[i] /= sum_freqs;
  // }

  return true;
}

bool
KmerProfile::indexPrefixes () {
  m_idx_pfx.clear();
  for (unsigned i=0; i<num_kmers; i++) {
    string pfx = m_vec_kmers[i].substr(0, kmer_length-1);
    if ( m_idx_pfx.find(pfx) == m_idx_pfx.end() ) {
      m_idx_pfx[pfx] = KmerProfile();
      m_idx_pfx[pfx].kmer_length = kmer_length;
    }
    m_idx_pfx[pfx].m_vec_kmers.push_back(m_vec_kmers[i]);
    m_idx_pfx[pfx].m_vec_weight.push_back(m_vec_weight[i]);
    m_idx_pfx[pfx].num_kmers++;
  }
  has_idx_pfx = true;
  return true;
}

} // namespace seqio