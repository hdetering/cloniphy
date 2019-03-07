#include "GenomeInstance.hpp"
#include "../stringio.hpp" // format()

using namespace std;
using boost::icl::interval_map;
using boost::icl::interval;
using stringio::format;

namespace seqio {

GenomeInstance::GenomeInstance () {}

GenomeInstance::GenomeInstance (
  const GenomeReference& g_ref
)
{
  for (auto const & kv : g_ref.chromosomes) {
    ChromosomeReference chr_ref = *(kv.second);
    // initial genome state is diploid -> generate two instances of each chromosome
    shared_ptr<ChromosomeInstance> sp_chr_inst1(new ChromosomeInstance(chr_ref, 'A'));
    this->vec_chr.push_back(sp_chr_inst1);
    this->vec_chr_len.push_back(sp_chr_inst1->length);
    shared_ptr<ChromosomeInstance> sp_chr_inst2(new ChromosomeInstance(chr_ref, 'B'));
    this->vec_chr.push_back(sp_chr_inst2);
    this->vec_chr_len.push_back(sp_chr_inst2->length);
    // sanity check: chromsome IDs should be unique
    assert (( format("Seqid '%s' already exists.", chr_ref.id.c_str()),
              this->map_id_chr.count(chr_ref.id)==0 ));
    // shared_ptr<ChromosomeInstance> sp_chr_inst1(up_chr_inst1);
    // shared_ptr<ChromosomeInstance> sp_chr_inst2(up_chr_inst2);
    this->map_id_chr[chr_ref.id] = { sp_chr_inst1, sp_chr_inst2 };
  }
}

GenomeInstance::GenomeInstance (
  const GenomeInstance& g_inst, 
  vector<seg_mod_t>& out_vec_seg_mod
)
{
  // copy chromosome lengths
  this->vec_chr_len = g_inst.vec_chr_len;

  // copy chromosome instances and index by name
  for (auto const & id_chr : g_inst.map_id_chr) {
    string id = id_chr.first;
    
    // update chromosome name index
    this->map_id_chr[id] = vector<shared_ptr<ChromosomeInstance>>();

    for (const shared_ptr<ChromosomeInstance> ci_old : id_chr.second) {
      ChromosomeInstance ci_new;
      ci_new.copy(ci_old, out_vec_seg_mod);
      this->vec_chr.push_back(make_shared<ChromosomeInstance>(ci_new));
      this->map_id_chr[id].push_back(make_shared<ChromosomeInstance>(ci_new));
    }
  }
}

void
GenomeInstance::addChromosome(
  shared_ptr<ChromosomeInstance> sp_chr,
  string id_chr)
{
  this->map_id_chr[id_chr].push_back(sp_chr);
  this->vec_chr.push_back(sp_chr);
  this->vec_chr_len.push_back(sp_chr->length);
}

void
GenomeInstance::deleteChromosome(
  shared_ptr<ChromosomeInstance> sp_chr,
  string id_chr)
{
  // remove chromosome from map
  this->map_id_chr[id_chr].erase(
    remove(
      this->map_id_chr[id_chr].begin(),
      this->map_id_chr[id_chr].end(),
      sp_chr),
    this->map_id_chr[id_chr].end());

  // rebuild indices
  this->vec_chr.clear();
  this->vec_chr_len.clear();
  for ( auto const & id_chr : this->map_id_chr ) {
    for ( auto const & sp_chr : id_chr.second ) {
      this->vec_chr.push_back(sp_chr);
      this->vec_chr_len.push_back(sp_chr->length);
    }
  }
}

void
GenomeInstance::duplicate (
  vector<seg_mod_t>& out_vec_seg_mod
)
{
  vector<tuple<string, shared_ptr<ChromosomeInstance>>> vec_tpl_id_ci;
  // loop over chromosome IDs in GenomeInstance
  for (auto const & kv : this->map_id_chr) {
    // copy all chromosomes under current chromosome ID
    for (auto const & ci_old : kv.second) {
      shared_ptr<ChromosomeInstance> ci_new(new ChromosomeInstance());
      ci_new->copy(ci_old, out_vec_seg_mod);
      vec_tpl_id_ci.push_back(make_tuple(kv.first, ci_new));
    }
  }
  // add newly created ChromosomeInstances to GenomeInstance
  for (auto tpl_id_chr : vec_tpl_id_ci) {
    string id_chr;
    shared_ptr<ChromosomeInstance> sp_chr;
    tie(id_chr, sp_chr) = tpl_id_chr;
    this->addChromosome(sp_chr, id_chr);
  }
}

vector<SegmentCopy>
GenomeInstance::getSegmentCopiesAt (
  string id_chr,
  TCoord ref_pos)
{
  vector<SegmentCopy> res_segments;
  // perform sanity checks
  assert( this->map_id_chr.count(id_chr) > 0 );

  // get SegmentCopies from all ChromosomeInstances
  for (auto const sp_chr : this->map_id_chr[id_chr]) {
    vector<SegmentCopy> chr_segs = sp_chr->getSegmentCopiesAt(ref_pos);
    res_segments.insert(res_segments.end(), chr_segs.begin(), chr_segs.end());
  }

  return res_segments;
}

void
GenomeInstance::getCopyNumberStates(
  map<unsigned, vector<shared_ptr<Locus>>>& map_cn_loci) {
  // for all chromosomes
  for ( auto const & id_chr : this->map_id_chr ) {
    string id = id_chr.first;
    interval_map<unsigned long, int> imap_reg_cn;
    // for all ChromosomeInstances
    for ( auto const & sp_chr : id_chr.second ) {
      // for all SegmentCopies
      for ( auto const & seg : sp_chr->lst_segments ) {
        auto i_reg = interval<unsigned long>::right_open(seg.ref_start, seg.ref_end);
        imap_reg_cn += make_pair(i_reg, 1);
      }
    }

    // collect genomic regions and CN states
    for (auto const & i_reg_cn : imap_reg_cn) {
      auto reg = i_reg_cn.first;
      int cn = i_reg_cn.second;
      if (map_cn_loci.count(cn) == 0)
        map_cn_loci[cn] = vector<shared_ptr<Locus>>();
      shared_ptr<Locus> sp_loc(new Locus(id, reg.lower(), reg.upper()));
      map_cn_loci[cn].push_back(sp_loc);
    }
  }
}

void
GenomeInstance::getCopyNumberStateByChr (
  map<string, interval_map<TCoord, double>>& map_chr_segments,
  const double scale
) const
{
  // infer copy number state segments
  // for each chromosome id, build an interval map
  for ( auto const & id_chr : map_id_chr ) {
    map_chr_segments[id_chr.first] = interval_map<TCoord, double>();
    // infer segment-wise copy number for each ChromosomeInstance
    for ( auto const chr : id_chr.second ) {
      // each SegmentCopy increases the CN state for the corresponding region
      for ( auto const & seg : chr->lst_segments ) {
        // NOTE: if this ever fails: switch start and end coordinates (or can interval_map deal with that?)
        assert( seg.ref_start < seg.ref_end );
        auto i = interval<TCoord>::right_open(seg.ref_start, seg.ref_end);
        // add interval with value 1 to interval map (copy number increases by 1)
        map_chr_segments[id_chr.first] += make_pair(i, scale);
      }
    }
  }
}

void
GenomeInstance::getCopyNumberStateByChr (
  map<string, interval_map<TCoord, AlleleSpecCopyNum>>& map_chr_alleles,
  const double scale
) const
{
  // infer copy number state segments
  // for each chromosome id, build an interval map
  for ( auto const & id_chr : map_id_chr ) {
    map_chr_alleles[id_chr.first] = interval_map<TCoord, AlleleSpecCopyNum>();
    // infer segment-wise copy number for each ChromosomeInstance
    for ( auto const chr : id_chr.second ) {
      // each SegmentCopy increases the CN state for the corresponding region
      for ( auto const & seg : chr->lst_segments ) {
        // NOTE: if this ever fails: switch start and end coordinates (or can interval_map deal with that?)
        assert( seg.ref_start < seg.ref_end );
        auto i = interval<TCoord>::right_open(seg.ref_start, seg.ref_end);
        // copy number increases for allele to which segment copy belongs
        AlleleSpecCopyNum allele_count;
        if ( seg.gl_allele == 'A' )
          allele_count.count_A = scale;
        else
          allele_count.count_B = scale;
        // add interval with value 1 to interval map (copy number increases by 1)
        map_chr_alleles[id_chr.first] += make_pair(i, allele_count);
      }
    }
  }
}


bool
GenomeInstance::indexSegmentCopies ()
{
  for (auto const id_chr : this->map_id_chr) {
    string id = id_chr.first;
    vector<shared_ptr<ChromosomeInstance>> vec_chr = id_chr.second;

    TSegMap imap_segments;
    for (const shared_ptr<ChromosomeInstance> sp_chr : vec_chr) {
      sp_chr->indexSegmentCopies(imap_segments);
    }
    this->map_chr_seg[id] = imap_segments;
  }

  return true;
}

ostream& operator<<(ostream& lhs, const GenomeInstance& gi) {
  lhs << "GenomeInstance" << endl;
  lhs << "--------------" << endl;
  for (auto const & id_chr : gi.map_id_chr) {
    unsigned num_ci = id_chr.second.size();
    lhs << "  Chromosome \"" << id_chr.first << "\" (" << num_ci << " instances):" << endl;
    for (auto const & ci : id_chr.second) {
      lhs << *ci;
    }
  }
  return lhs;
}

} // namespace seqio