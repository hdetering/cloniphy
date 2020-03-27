#include "BulkSample.hpp"

using namespace std;
using boost::icl::interval;
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
  const string chr,
  const TCoord pos,
  const double cvg_per_bp
) const
{
  double cvg_exp = -1.0;

  // get total copy number state for variant locus
  double seg_cn  = -1.0; 
  TCoord seg_len = 0;
  this->getTotalCopyNumberAt(chr, pos, seg_cn, seg_len);
  
  if (seg_cn < 0) { // this is bad!
    fprintf(stderr, "[ERROR] (BulkSample::getExpectedCoverageAt) Failed to get copy number state for '%s:%lu'\n", 
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
  double& out_cn_tot,
  TCoord& out_seg_len
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

bool
BulkSample::initAlleleCounts (
  const map<string, double> map_clone_ccf,
  const vario::VariantStore& var_store,
  const map<string, map<string, seqio::TSegMap>> map_clone_chr_seg
  //vario::TMapChrPosVaf& out_map_chr_pos_vaf
)
{
  // requires genome instances to be initialized
  assert(has_clone_genomes);

  this->m_map_clone_snv_vac.clear();
  this->m_map_snv_vaf.clear();

  // loop over clone genomes
  for (auto clone_chr_seg : map_clone_chr_seg) {

    string id_clone = clone_chr_seg.first;
    // segment copy map for clone, indexed by chromosome
    map<string, seqio::TSegMap> map_chr_seg = clone_chr_seg.second;
    // clone cellular prevalence in this BulkSample
    double clone_weight = map_clone_ccf.at(id_clone);

    // initialize map for clone
    m_map_clone_snv_vac[id_clone] = map<int, vario::VariantAlleleCount>();

    //--------------------------------------------------------------------------
    // populate allele counts for SNVs
    //--------------------------------------------------------------------------
    for (auto const & kv : var_store.map_id_snv) {
      int id_var = kv.first;
      Variant var = kv.second;
      
      TCoord pos_var = var.pos;
      string id_chr = var.chr;

      // if (out_map_chr_pos_vaf.count(id_chr) == 0) {
      //   out_map_chr_pos_vaf[id_chr] = TMapPosVaf();
      // }
      
      // get segment copies overlapping SNV position
      seqio::TSegMap imap_seg_chr = map_clone_chr_seg.at(id_clone).at(id_chr);
      seqio::TSegSet iset_read;
      iset_read.add(interval<TCoord>::right_open(pos_var, pos_var+1));
      seqio::TSegMap imap_seg_ovlp = imap_seg_chr & iset_read;

      short num_tot = 0;
      short num_alt = 0;

      for (auto const & itvl_iset : imap_seg_ovlp) {
        for (const SegmentCopy & segment : itvl_iset.second) {
          num_tot++; // increase total allele count
          // check if current segment copy carries current SNV
          auto it_seg_vars = var_store.map_seg_vars.find(segment.id);
          if (it_seg_vars == var_store.map_seg_vars.end()) 
            continue;
          vector<int> vec_seg_vars = it_seg_vars->second;

          // increase alternative allele count if variant associated to segment copy
          if(find(vec_seg_vars.begin(), vec_seg_vars.end(), id_var) != vec_seg_vars.end())
            num_alt++;
        }
      }

      vario::VariantAlleleCount vac;
      vac.num_tot = num_tot;
      vac.num_alt = num_alt;
      m_map_clone_snv_vac[id_clone][id_var] = vac;
    }
  }

  // Initialize expected bulk variant allele frequencies (VAFs) for somatic SNVs.
  // Formula:
  //   VAF_i = \sum_k ( CP_k * A_k,i / N_k,i )
  // With
  //   CCF_k := Cancer cell fraction (frequency) of clone k within the bulk sample
  //   A_k,i := Alternative alleles of clone k at variant i
  //   N_k,i := Total alleles of clone k at variant i

  for (auto const & kv : var_store.map_id_snv) {
    int id_snv = kv.first;
    double vaf = 0.0;

    for (auto const & clone_ccf : map_clone_ccf) {
      string id_clone = clone_ccf.first;
      double ccf = clone_ccf.second;
      vario::VariantAlleleCount vac = this->m_map_clone_snv_vac[id_clone][id_snv];

      vaf += ccf * vac.num_alt / vac.num_tot;
    }

    this->m_map_snv_vaf[id_snv] = vaf;
  }
}

} // namespace bamio