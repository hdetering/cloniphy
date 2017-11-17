#include "BulkSample.hpp"

using namespace std;
using boost::icl::interval_map;
using seqio::TCoord;
using seqio::AlleleSpecCopyNum;

namespace bamio {


BulkSample::BulkSample() {}

BulkSample::BulkSample(
  const string id_sample,
  const map<string, double> clone_weight
)
: id(id_sample),
  m_clone_weight(clone_weight),
  genome_len_abs(0)
{}

vector<string>
BulkSample::getCloneLabels () const
{
  vector<string> vec_clone_lbl;

  for (auto lbl_w : this->m_clone_weight) {
    vec_clone_lbl.push_back(lbl_w.first);
  }

  return vec_clone_lbl;
}

double
BulkSample::getExpectedCoverageAt (
  string chr,
  TCoord pos,
  double cvg_per_bp
)
{
  double cvg_exp = -1.0;

  // get total copy number state for variant locus
  double seg_cn  = -1.0; 
  TCoord seg_len = 0;
  this->getTotalCopyNumberAt(chr, pos, seg_cn, seg_len);
  
  if (seg_cn < 0) { // this is bad!
    fprintf(stderr, "[ERROR] Failed to get copy number state for '%s:%lu'", 
      chr.c_str(), 
      pos);
  } 
  else {
    cvg_exp = cvg_per_bp * seg_cn * seg_len;
  }

  return cvg_exp;
}

bool
BulkSample::getTotalCopyNumberAt (
  const string chr,
  const TCoord pos,
  double out_cn_tot,
  TCoord out_seg_len
) const
{
  typedef interval_map<TCoord, AlleleSpecCopyNum> TCnMap;

  // sanity check: chromosome id present? 
  if ( this->m_chr_cn.count(chr) == 1 ) {
    // find interval overlapping position
    TCnMap::const_iterator it = this->m_chr_cn.at(chr).find(pos);
    // retrieve interval that contains given locus
    auto seg = it->first;
    out_seg_len = seg.upper() - seg.lower();
    // retrieve allele-specific (maternal/paternal) copy number state
    AlleleSpecCopyNum ascn = it->second;
    // total copy number: sum of maternal and paternal allele counts
    out_cn_tot = ascn.count_A + ascn.count_B;
  }
  else { // this is to be treated as an error
    return false;
  }

  return true;
}


} // namespace bamio