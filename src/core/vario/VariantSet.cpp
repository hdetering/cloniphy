#include "VariantSet.hpp"

using namespace std;

namespace vario {

  /* VariantSet
 *------------*/

VariantSet::VariantSet() {}
VariantSet::~VariantSet() {}

VariantSet::VariantSet(vector<Variant> variants) {
  this->vec_variants = variants;
  this->indexVariants();
  this->calculateSumstats();
}

VariantSet& VariantSet::operator+=(const VariantSet& rhs) {
  this->vec_variants.insert(
    this->vec_variants.end(),
    rhs.vec_variants.begin(),
    rhs.vec_variants.end()
  );
  this->indexVariants();
  this->calculateSumstats();
  return *this; // return the result by reference
}

long VariantSet::indexVariants() {
  long num_variants;
  this->map_chr2pos2var.clear();

  for (Variant var : this->vec_variants) {
    if (this->map_chr2pos2var.find(var.chr) == this->map_chr2pos2var.end())
      this->map_chr2pos2var[var.chr] = map<unsigned long, vector<Variant>>();
    //if (this->map_chr2pos2var[var.chr].find(var.pos) == this->map_chr2pos2var[var.chr].end())
    //  this->map_chr2pos2var[var.chr][var.pos] = vector<Variant>>();
    this->map_chr2pos2var[var.chr][var.pos].push_back(var);
    num_variants++;
  }

  this->num_variants = num_variants;
  return num_variants;
}

long VariantSet::calculateSumstats() {
  string nucs = "ACGTacgt";
  long num_variants = 0;
  long num_substitutions = 0;
  long mat_subs[4][4] = {
    { 0, 0, 0, 0},
    { 0, 0, 0, 0},
    { 0, 0, 0, 0},
    { 0, 0, 0, 0}
  };

  for (Variant v : vec_variants) {
    num_variants++;
    string ref = v.alleles[0];
    if (ref.length() == 1) { // make sure we're not dealing with an InDel
      char ref_nuc = ref[0];
      for (unsigned i=1; i<v.alleles.size(); ++i) {
        string alt = v.alleles[i];
        if (alt.length() == 1) { // make sure we're not dealing with an InDel
          char alt_nuc = alt[0];
          if (nucs.find(alt_nuc) != string::npos) {
            num_substitutions++;
            short ref_idx = seqio::nuc2idx(ref_nuc);
            short alt_idx = seqio::nuc2idx(alt_nuc);
            mat_subs[ref_idx][alt_idx] += 1;
          }
        }
      }
    }
  }

  // normalize substitution counts
  for (auto i=0; i<4; ++i) {
    for (auto j=0; j<4; ++j) {
      mat_freqs[i][j] = (double)mat_subs[i][j] / num_substitutions;
    }
  }

  this->num_variants = num_variants;
  return num_substitutions;
}


} // namespace vario