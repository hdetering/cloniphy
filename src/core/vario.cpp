#include "vario.hpp"
#include <algorithm>
#include <boost/container/flat_set.hpp>
#include <boost/format.hpp>
#include <ctime>
#include <fstream>
#include <set>
#include <stdio.h>

using namespace std;
using boost::format;
using boost::uuids::uuid;
using seqio::ChromosomeInstance;
using seqio::Locus;
using seqio::Nuc;

namespace vario {

Variant::Variant(std::string id, std::string chr, unsigned long pos)
 : id(id), chr(chr), pos(pos) {}

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

/** Index variants by chromosome and position. */
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

Mutation::Mutation() {
  this->id = 0;
  this->is_snv = false;
  this->is_cnv = false;

  //TODO: deprecated
  // this->relPos = 0.0;
  // this->copy = 0;
}

/*
Mutation::Mutation(char ref, char alt) {
  short ref_nuc = seqio::charToNuc(ref);
  short alt_nuc = seqio::charToNuc(alt);
  this->offset =  ((alt_nuc-ref_nuc) % 4); // TODO: this could be more generic (get rid of the hard-coded 4)
}*/

// bool Mutation::operator< (const Mutation &other) const {
//   return (relPos + copy) < (other.relPos + other.copy);
// }

// vector<Mutation> Mutation::sortByPosition(const vector<Mutation> &mutations) {
//   vector<Mutation> mutationsCopy = mutations;
//   sort(mutationsCopy.begin(), mutationsCopy.end());
//   return mutationsCopy;
// }

CopyNumberVariant::CopyNumberVariant ()
: id(0),
  ref_chr(""),
  ref_pos_begin(0),
  ref_pos_end(0),
  is_deletion(false),
  is_wgd(false)
{}

Variant::Variant ()
: id(""),
  chr(""),
  pos(0),
  alleles(0),
  idx_mutation(0),
  rel_pos(0.0),
  is_somatic(false)
{}
Variant::~Variant () {}

bool Variant::operator< (const Variant &other) const {
  return rel_pos < other.rel_pos;
}

vector<Variant> Variant::sortByPosition(const vector<Variant> &variants) {
  vector<Variant> variantsCopy = variants;
  sort(variantsCopy.begin(), variantsCopy.end());
  return variantsCopy;
}

vector<Variant> Variant::sortByPositionLex(const vector<Variant> &variants) {
  vector<Variant> variantsCopy = variants;
  sort(variantsCopy.begin(), variantsCopy.end(),
      [](const Variant &a, const Variant &b) -> bool {
        if (a.chr < b.chr) return true;
        if (a.chr == b.chr) return a.pos < b.pos;
        return false;
      });
  return variantsCopy;
}

vector<Variant> Variant::sortByPositionPoly(const vector<Variant> &variants) {
  vector<Variant> variantsCopy = variants;
  sort(variantsCopy.begin(), variantsCopy.end(),
      [](const Variant &a, const Variant &b) -> bool {
        return (a.rel_pos + a.chr_copy) < (b.rel_pos + b.chr_copy);
      });
  return variantsCopy;
}

vector<Variant> Variant::sortByPositionRef(const vector<Variant> &variants) {
  vector<Variant> variantsCopy = variants;
  sort(variantsCopy.begin(), variantsCopy.end(),
      [](const Variant &a, const Variant &b) -> bool {
        return a.pos < b.pos;
      });
  return variantsCopy;
}

bool Variant::isSnv() {
  for (vector<string>::iterator allele=this->alleles.begin(); allele!=this->alleles.end(); ++allele) {
    if ((*allele).size()>1) { return false; }
  }
  return true;
}

/*------------------------------------*/
/*           Utility methods          */
/*------------------------------------*/

unsigned assignMutationType(
  std::vector<Mutation>& vec_mutations,
  const double ratio_cnv,
  RandomNumberGenerator<>& rng)
{
  function<double()> random_double = rng.getRandomFunctionDouble(0.0, 1.0);
  unsigned num_mut=0;
  unsigned num_cnv=0;

  for (auto &m : vec_mutations) {
    double r_dbl = random_double();
//fprintf(stderr, "<Mutation(id=%u;relPos=%f,copy=%d)>\n", i, rel_pos, copy);
    m.id = num_mut++;
    if (r_dbl < ratio_cnv) {
      m.is_snv = false;
      m.is_cnv = true;
      num_cnv++;
    } else {
      m.is_snv = true;
      m.is_cnv = false;
    }
  }

  return num_cnv;
}

/** Read mutation map (clone x mutation) from a CSV file. */
int readMutMapFromCSV(
  map<string, vector<bool>> &mm,
  const string &filename)
{
  unsigned n_rows_mm = 0;
  vector<vector<string>> mtx_csv;
  n_rows_mm = stringio::readCSV(mtx_csv, filename);
  // convert CSV table to mutation map
  mm.clear();
  for (auto row : mtx_csv) {
    string c_lbl = row[0];
    unsigned n_col = row.size();

    vector<bool> vec_row(n_col-1, false);
    for (auto i=1; i<n_col; i++) {
      vec_row[i-1] = (row[i] == "0" ? false : true);
    }

    mm[c_lbl] = vec_row;
  }
}

void readVcf(
  string fn_vcf,
  VariantSet& variants,
  map<string, vector<Genotype> > &gtMatrix)
{
  ifstream f_vcf;
  f_vcf.open(fn_vcf.c_str(), ios::in);
  readVcf(f_vcf, variants, gtMatrix);
  f_vcf.close();
}

void readVcf(
  istream &input,
  VariantSet& variants,
  map<string, vector<Genotype> > &gtMatrix)
{
  unsigned num_samples = 0;
  string line = "";
  string header_line;
  // consume header lines
  while (stringio::safeGetline(input, line) && line[0]=='#') {
    header_line = line;
  }
fprintf(stderr, "VCF header: %s\n", header_line.c_str());
  // parse header
  vector<string> hdr_cols = stringio::split(header_line, '\t');
  vector<string> lbl_samples(hdr_cols.begin()+9, hdr_cols.end()); // samples start at column 10
  num_samples = lbl_samples.size();
  gtMatrix.clear();
  for (auto lbl : lbl_samples) {
    gtMatrix[lbl] = vector<Genotype>();
  }

  // parse variants
  unsigned var_idx = 0;
  while(input.good()) {
    vector<string> var_cols = stringio::split(line, '\t');
    if (var_cols.size() != num_samples+9) {
      fprintf(stderr, "[ERROR] number of columns does not match header (variant '%s').\n", var_cols[2].c_str());
      stringio::safeGetline(input, line);
      continue;
    }
    Variant var;
    var.chr = var_cols[0];
    var.pos = atol(var_cols[1].c_str())-1; // NOTE: internal positions should always zero-based!
    var.id  = var_cols[2];
    var.alleles.push_back(var_cols[3]);
    vector<string> alt = stringio::split(var_cols[4], ',');
    var.alleles.insert(var.alleles.end(), alt.begin(), alt.end());
    var.reg_copy = 0; // TODO: can we do better than this?
    var.chr_copy = 0; // TODO: can we do better than this?
    // TODO: at this point only SNVs are supported
    if (!var.isSnv()) {
      stringio::safeGetline(input, line);
      continue;
    }

    // for single-sample VCFs, set variant phase
    if (num_samples == 1) {
      if (atoi(&var_cols[9][0]) > 0) {
        var.chr_copy = 0;
      } else if (atoi(&var_cols[9][2]) > 0) {
        var.chr_copy = 1;
      }
    }

    variants.vec_variants.push_back(var);
    var_idx++;
//fprintf(stderr, "read from VCF: <Variant(idx=%u,id=%s,pos=%lu,num_alleles=%lu)> ", var_idx, var->id.c_str(), var->pos, var->alt.size()+1);
    for (unsigned i=0; i<num_samples; ++i) {
      string genotype = var_cols[9+i];
      short a_allele = atoi(&genotype[0]);
      short b_allele = atoi(&genotype[2]); // TODO: this will fail for >9 alleles.
//fprintf(stderr, "Genotype for sample %u: %d, %d\n", i, a_allele, b_allele);
      Genotype gt = { var.id, a_allele, b_allele };
      gtMatrix[lbl_samples[i]].push_back(gt);
    }
    stringio::safeGetline(input, line);
  }
  variants.num_variants = var_idx;
  // index Variants by chromosome and position
  variants.indexVariants();
  // calculate substitution freqs etc.
  variants.calculateSumstats();
}

void writeVcf(
  const vector<shared_ptr<SeqRecord>>& seqs,
  const vector<Variant>& vars,
  const vector<int>& id_samples,
  const vector<string>& labels,
  const vector<vector<bool> >& mutMatrix,
  ostream& out)
{
  unsigned num_samples = id_samples.size();
  string var_qual = "."; // QUAL column
  string var_info = (format("NS=%d") % num_samples).str(); // INFO column
  string var_fmt  = "GT:GQ"; // FORMAT column
  short  gt_qual  = 60; // genotype quality (global)

  // write header
  out << "##fileformat=VCFv4.1" << endl;
  time_t timer = time(NULL);
  tm* t = localtime(&timer);
  out << format("##fileDate=%d-%d-%d") % (1900+t->tm_year) % t->tm_mon % t->tm_mday << endl;
  out << "##source=CloniPhy v0.01" << endl;
  //out << "##reference=" << std::endl; # TODO: include ref filename
  vector<string> vec_ref_ids;
  for (auto rec : seqs) {
    if (find(vec_ref_ids.begin(), vec_ref_ids.end(), rec->id_ref) == vec_ref_ids.end()) {
      vec_ref_ids.push_back(rec->id_ref);
      out << format("##contig=<ID=%d, length=%u>") % rec->id_ref % rec->seq.size() << endl;
    }
  }
  out << "##phasing=complete" << endl;
  out << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" << endl;
  out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
  out << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << endl;
  out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  //out << "\thealthy";
  for (auto lbl : labels)
    out << '\t' << lbl;
  out << endl;

  // write variants
  vector<Variant> sorted_vars = Variant::sortByPositionLex(vars);
  for (Variant var : sorted_vars) {
    string ref = var.alleles[0];
    string alt = var.alleles[1];
    out << format("%s\t%d\t%d\t%s\t%s\t%s\tPASS\t%d\t%s")
                  % (var.chr) % (var.pos+1) % (var.id)
                  % ref % alt % var_qual % var_info % var_fmt;
    for (auto sid : id_samples) {
      string genotype = "";
      if (mutMatrix[sid][var.idx_mutation] == true)
        genotype = (var.chr_copy==0 ? "1|0" : "0|1");
      else
        genotype = "0|0";
      out << format("\t%s:%d") % genotype % gt_qual;
    }
    out << endl;
  }
}

/** Generate VCF output from a reference genome and a set of variants.
    (single sample) */
void writeVcf(
  const vector<shared_ptr<SeqRecord>>& seqs,
  const vector<Variant>& vars,
  const string label,
  const string filename)
{
  vector<int> vec_ids(1, 0);
  vector<string> vec_labels(1, label);
  vector<vector<bool>> mtx_mut(1, vector<bool>(vars.size(), true));

  ofstream f_out;
  f_out.open(filename);
  writeVcf(seqs, vars, vec_ids, vec_labels, mtx_mut, f_out);
  f_out.close();
}

long getRandomMutPos(
  function<int()> r_bucket,
  map<string, vector<long>> map_bucket_pos,
  short offset)
{
  return 0;
}

/** Generate variant loci in a given genome based on evolutionary model.
    Nucleotide substitution probabilities guide selection of loci. */
vector<Variant> generateGermlineVariants(
  const int num_variants,
  const GenomeReference& genome,
  GermlineSubstitutionModel& model,
  RandomNumberGenerator<>& rng,
  const bool inf_sites)
{
  vector<Variant> variants = vector<Variant>(num_variants);
  boost::container::flat_set<int> var_pos; // keep track of variant positions
  function<double()> random_float = rng.getRandomFunctionDouble(0.0, 1.0);
  function<short()> random_copy = rng.getRandomFunctionInt(short(0), short(genome.ploidy-1));
  random_selector<> selector(rng.generator); // used to pick random vector indices

  // determine base mutation probs from model (marginal sums)
  vector<double> p_i(4, 0);
  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      p_i[i] += model.Qij[i][j];
    }
  }
  function<int()> random_nuc_idx = rng.getRandomIndexWeighted(p_i);

  unsigned long genome_len = genome.length; // haploid genome length
  for (int i=0; i<num_variants; ++i) {
    // pick random nucleotide bucket
    int idx_bucket = random_nuc_idx();
    // pick random position
    long nuc_pos = selector(genome.nuc_pos[idx_bucket]);
    if (inf_sites) {
      while (binary_search(var_pos.begin(), var_pos.end(), nuc_pos)) {
// TODO: check verbosity setting
fprintf(stderr, "[INFO] Infinite sites assumption: locus %ld has been mutated before, picking another one...\n", nuc_pos);
        nuc_pos = selector(genome.nuc_pos[idx_bucket]);
      }
      var_pos.insert(nuc_pos);
    }
    Locus loc = genome.getLocusByGlobalPos(nuc_pos);
    // pick new nucleotide
    short nuc_alt = evolution::MutateSite(idx_bucket, random_float, model);
    Variant var;
    var.id = str(format("g%d") % i);
    var.chr = loc.id_ref;
    var.chr_copy = random_copy();
    var.rel_pos = double(nuc_pos-(var.chr_copy*genome_len))/genome_len;
    var.reg_copy = 0; // TODO: when implementing CNVs, use this property to indicate affected copy
    var.pos = loc.start;
    var.alleles.push_back(string(1, seqio::idx2nuc(idx_bucket)));
    var.alleles.push_back(string(1, seqio::idx2nuc(nuc_alt)));
    var.idx_mutation = i;
    variants[i] = var;
  }

  return variants;
}

/** Generate variant loci in a given genome based on somatic mutation model. */
void VariantStore::generateSomaticVariants(
  const vector<Mutation>& vec_mutations,
  const GenomeReference& genome,
  SomaticSubstitutionModel& model_snv,
  SomaticCnvModel& model_cnv,
  RandomNumberGenerator<>& rng,
  const bool inf_sites)
{
  // SNV events
  //------------
  vector<Variant> variants;
  // keep track of variant positions (ISM)
  boost::container::flat_set<int> var_pos;
  // random function, returns substitution index
  function<int()> r_idx_sub = rng.getRandomIndexWeighted(model_snv.m_weight);
  function<short()> random_copy = rng.getRandomFunctionInt(short(0), short(genome.ploidy-1));
  random_selector<> selector(rng.generator); // used to pick random vector indices
  unsigned long genome_len = genome.length; // haploid genome length
  //------------

  // CNV events
  //------------
  // random function, selects CNV event type
  // random function, selects CNV event length (relative to CHR len)
  // chromosome lengths (used to select affected chromosome by size)
  //------------

  for (auto m : vec_mutations) {
    if (m.is_snv) {
      // pick random ref->alt substitution
      int i_sub = r_idx_sub();
      string ref_site = model_snv.m_site[i_sub];
      string alt_nuc = model_snv.m_alt[i_sub];
      string ref_nuc = ref_site.substr(1, 1);
      // pick random position (+1 b/c second nucleotide in 3-mer is mutated)
      long nuc_pos = selector(genome.map_3mer_pos.at(ref_site)) + 1;
      if (inf_sites) {
        while (binary_search(var_pos.begin(), var_pos.end(), nuc_pos)) {
          // TODO: check verbosity setting
          fprintf(stderr, "[INFO] Infinite sites model: locus %ld has been mutated before, picking another one...\n", nuc_pos);
          nuc_pos = selector(genome.map_3mer_pos.at(ref_site)) + 1;
        }
        var_pos.insert(nuc_pos);
      }
      Locus loc = genome.getLocusByGlobalPos(nuc_pos);
      // TODO: identify available segment copies in GenomeInstance, choose one

      // init new Variant
      Variant var;
      var.id = str(format("s%d") % m.id);
      var.chr = loc.id_ref;
      var.chr_copy = random_copy();
      var.rel_pos = double(nuc_pos-(var.chr_copy*genome_len))/genome_len;
      // TODO: assign SegmentCopy
      //var.reg_copy = 0;
      var.pos = loc.start;
      var.alleles.push_back(ref_nuc);
      var.alleles.push_back(alt_nuc);
      var.idx_mutation = m.id;
      var.is_somatic = true;
      this->map_id_snv[m.id] = var;
      // TODO: add variant to map_seg_vars
    }
    else {
      // TODO: simulate CNV events
      // pick event type
      // if WGD:

      //
    }
  }

  //return variants;
}

void VariantStore::applyMutation(
  Mutation mut,
  GenomeInstance& genome,
  RandomNumberGenerator<>& rng)
{
  // perform some sanity checks
  assert( mut.is_snv != mut.is_cnv );
  assert( !mut.is_snv || this->map_id_snv.count(mut.id)>0 );
  assert( !mut.is_cnv || this->map_id_cnv.count(mut.id)>0 );

  random_selector<> selector(rng.generator); // used to pick random SegmentCopy

  if ( mut.is_snv ) { // SNV mutation
    Variant snv = this->map_id_snv[mut.id];
    // get available SegmentCopies
    vector<SegmentCopy> seg_targets = genome.getSegmentCopiesAt(snv.chr, snv.pos);
    SegmentCopy sc = selector(seg_targets);
    // initialize or append to Variant vector of SegmentCopy
    if ( map_seg_vars.count(sc.id)==0 ) {
      map_seg_vars[sc.id] = { mut.id };
    } else {
      map_seg_vars[sc.id].push_back(mut.id);
    }
  }
  else { // CNV mutation
    // TODO: apply copy number variant
    CopyNumberVariant cnv = this->map_id_cnv[mut.id];
    // track modifications to SegmentCopies
    vector<seqio::seg_mod_t> vec_seg_mod;

    if (cnv.is_wgd) { // Whole Genome Duplication
      genome.duplicate(vec_seg_mod);
      this->transferMutations(vec_seg_mod);
    }
    else { // Not a WGD, chromosome region is affected
      // pick a ChromsomeInstance randomly (weighted by chromosome lengths)
      string id_chr = cnv.ref_chr;
      assert( genome.map_id_chr.find(id_chr) != genome.map_id_chr.end() );
      vector<unsigned long> chr_len;
      for (auto const & chr_inst : genome.map_id_chr[id_chr]) {
        chr_len.push_back(chr_inst->length);
      }
      function<int()> r_idx_chr = rng.getRandomIndexWeighted(chr_len);
      int idx_chr = r_idx_chr();
      shared_ptr<ChromosomeInstance> sp_chr = genome.map_id_chr[id_chr][idx_chr];

      if (cnv.is_deletion) { // delete a genomic region
        vec_seg_mod = sp_chr->deleteRegion(cnv.start_rel, cnv.len_rel, cnv.is_forward, cnv.is_telomeric);
        this->transferMutations(vec_seg_mod);
      } else { // amplify a genomic region
        vec_seg_mod = sp_chr->amplifyRegion(cnv.start_rel, cnv.len_rel, cnv.is_forward, cnv.is_telomeric);
        this->transferMutations(vec_seg_mod);
      }
    }
  }
}

void VariantStore::transferMutations(vector<seqio::seg_mod_t> vec_seg_mod) {
  uuid seg_new_id;
  uuid seg_old_id;
  unsigned long seg_old_start;
  unsigned long seg_old_end;

  // loop over SegmentCopy modifications
  for (auto const & tpl_seg_mod : vec_seg_mod) {
    // each modification carries information about:
    // - newly created SegmentCopy
    // - existing SegmentCopy it was copied from
    // - interval (start, end) the new SegmentCopy originates from
    tie(seg_new_id, seg_old_id, seg_old_start, seg_old_end) = tpl_seg_mod;

    // transfer variants associated with the copied region within old SegmentCopy
    vector<unsigned> seg_new_vars;
    auto it_vars = this->map_seg_vars.find(seg_old_id);
    if (it_vars == this->map_seg_vars.end())
      continue;
    for (auto const & id_snv : it_vars->second) {
      Variant snv = this->map_id_snv[id_snv];
      if (snv.pos >= seg_old_start && snv.pos < seg_old_end) {
        seg_new_vars.push_back(id_snv);
      }
    }
    if (seg_new_vars.size() > 0) {
      this->map_seg_vars[seg_new_id] = seg_new_vars;
    }
  }
}


vector<Variant> VariantStore::getSnvVector() {
  vector<Variant> variants;
  for (auto const & kv : this->map_id_snv)
    variants.push_back(kv.second);

  return variants;
}

/** Generate variant loci in a given genome based on evolutionary model.
    Nucleotide substitution probabilities guide selection of loci. */
vector<Variant> generateVariantsRandomPos(
  const int num_variants,
  const GenomeReference& genome,
  GermlineSubstitutionModel& model,
  RandomNumberGenerator<>& rng,
  const bool inf_sites)
{
  vector<Variant> variants = vector<Variant>(num_variants);
  boost::container::flat_set<int> var_pos; // keep track of variant positions
  function<double()> random_float = rng.getRandomFunctionDouble(0.0, 1.0);
  function<long()> random_pos = rng.getRandomFunctionInt<long>(0, genome.length);
  function<short()> random_copy = rng.getRandomFunctionInt(short(0), short(genome.ploidy-1));
  random_selector<> selector(rng.generator); // used to pick random vector indices

  for (int i=0; i<num_variants; ++i) {
    // pick random position
    long nuc_pos = random_pos();
    if (inf_sites) {
      while (binary_search(var_pos.begin(), var_pos.end(), nuc_pos)) {
fprintf(stderr, "locus %ld has beend mutated before, picking another one...\n", nuc_pos);
        nuc_pos = random_pos();
      }
      var_pos.insert(nuc_pos);
    }
    Locus loc = genome.getLocusByGlobalPos(nuc_pos);
    string id_chr = genome.records[loc.idx_record]->id;
    short ref_nuc = seqio::nuc2idx(genome.records[loc.idx_record]->seq[loc.start]);
    // pick new nucleotide
    short nuc_alt = evolution::MutateSite(ref_nuc, random_float, model);
    Variant var;
    var.id = to_string(i);
    var.chr = id_chr;
    var.chr_copy = random_copy();
    var.reg_copy = 0; // TODO: when implementing CNVs, use this property to indicate affected copy
    var.pos = nuc_pos;
    var.rel_pos = double(nuc_pos)/genome.length;
    var.alleles.push_back(string(1, seqio::idx2nuc(ref_nuc)));
    var.alleles.push_back(string(1, seqio::idx2nuc(nuc_alt)));
    variants[i] = var;
  }

  return variants;
}

void applyVariants(
  GenomeReference &genome,
  const vector<Variant> &variants)
{
  unsigned num_sequences = genome.num_records;
  // generate lookup table for sequences
  map<string,vector<unsigned>> chr2seq;
  for (unsigned i=0; i<genome.records.size(); ++i) {
    string id_ref = genome.records[i]->id_ref;
    if (chr2seq.find(id_ref) == chr2seq.end()) {
      chr2seq[id_ref] = vector<unsigned>();
    }
    chr2seq[id_ref].push_back(i);
  }
fprintf(stderr, "applying %lu variants...\n", variants.size());
  // modify genome according to variant genotypes
  for (Variant var : variants) {
    auto it_chr_idx = chr2seq.find(var.chr);
    // make sure ref seq exists in genome
    if (it_chr_idx == chr2seq.end()) {
      fprintf(stderr, "[WARN] sequence '%s' not contained in reference genome\n", var.chr.c_str());
    }
    else {
      vector<unsigned> chr_idx = it_chr_idx->second;
      // make sure mutated chr copy exists in genome
      if (chr_idx.size() < var.chr_copy+1) {
        fprintf(stderr, "[WARN] genome does not contain %d copies of CHR '%s'\n", var.chr_copy+1, var.chr.c_str());
      }
      else {
        // apply variant to sequence
        unsigned cidx = chr_idx[var.chr_copy];
        genome.records[cidx]->seq[var.pos] = var.alleles[1][0]; // TODO: at the moment only SNVs are supported ("[0]" extracts the first character from the allel)
      }
    }
  }
}


void applyVariants(
  GenomeReference &genome,
  const vector<Variant> &variants,
  const vector<Genotype> &genotypes)
{
  unsigned num_sequences = genome.num_records;
  // generate lookup table for sequences
  map<string,unsigned> chr2seq;
  for (unsigned i=0; i<num_sequences; ++i) {
    chr2seq[genome.records[i]->id_ref] = i;
  }
fprintf(stderr, "variants: %lu\n", variants.size());
fprintf(stderr, "genotypes: %lu\n", genotypes.size());
  // modify genome according to variant genotypes
  for (unsigned i=0; i<variants.size(); ++i) {
    Variant var = variants[i];
//fprintf(stderr, "variant %u\n", i);
    Genotype gt = genotypes[i];
    map<string,unsigned>::iterator it_chr_idx = chr2seq.find(var.chr);
    // make sure ref seq exists in genome
    if (it_chr_idx != chr2seq.end()) {
      unsigned chr_idx = it_chr_idx->second;
      // only apply variant if any allele is non-reference
      if (gt.maternal>0) {
        genome.records[chr_idx]->seq[var.pos-1] = var.alleles[gt.maternal][0]; // TODO: at the moment only SNVs are supported ("[0]" extracts the first character from the allel)
      }
      if (gt.paternal>0) {
//fprintf(stderr, "\t%lu paternal: '%s' -> '%s'\n", var.pos, var.alleles[0].c_str(), var.alleles[gt.paternal].c_str());
        genome.records[num_sequences+chr_idx]->seq[var.pos-1] = var.alleles[gt.paternal][0];
      }
    }
    else {
      fprintf(stderr, "[WARN] sequence '%s' not contained in reference genome\n", var.chr.c_str());
      return;
    }
  }
}

// void applyVariantsStream(
//   const Genome &g,
//   const vector<Mutation> &mutations,
//   const vector<Variant> &variants,
//   ostream &outstream,
//   short len_line)
// {
//   // sort mutations according to (diploid) genomic position
//   vector<Mutation> vec_mut_sorted = Mutation::sortByPosition(mutations);
//   unsigned idx_mut = 0;
//   unsigned idx_chr = 0;
//   unsigned idx_nuc = 0;
//   short idx_char = 0;
//   bool chr_mutated = false;
//   Variant next_var = variants[vec_mut_sorted[idx_mut].id];
//   for (idx_chr=0; idx_chr<g.records.size(); ++idx_chr) {
//     chr_mutated = (g.records[idx_chr].id_ref == next_var.chr); // is this chromosome mutated?
//     chr_mutated = chr_mutated && (g.records[idx_chr].chr_copy == vec_mut_sorted[idx_mut].copy); // is this sequence the right copy?
//     // print header
//     outstream << ">" << g.records[idx_chr].id << endl;
//     idx_nuc = 0;
//     for (string::const_iterator nuc=g.records[idx_chr].seq.begin();
//          nuc!=g.records[idx_chr].seq.end(); ++nuc) {
//       if (chr_mutated && idx_nuc++ == next_var.pos) {
//         outstream << next_var.alleles[1]; // print variant nucleotide
//         if (++idx_mut<vec_mut_sorted.size())
//           next_var = variants[vec_mut_sorted[idx_mut].id];
//         else
//           chr_mutated = false;
//       }
//       else // print reference nucleotide
//         outstream << *nuc;
//       if (++idx_char == len_line) { // enforce fixed line width
//         outstream << endl;
//         idx_char = 0;
//       }
//     }
//     if (idx_char != 0) { // end of sequence
//       outstream << endl;
//       idx_char = 0;
//     }
//   }
// }

} /* namespace vario */
