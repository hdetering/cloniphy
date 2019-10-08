#include "vario.hpp"
//#include "seqio/ChromosomeInstance.hpp"
#include <algorithm>
#include <boost/container/flat_set.hpp>
#include <ctime>
#include <fstream>
#include <set>
#include <cstdio>

using namespace std;
using boost::uuids::uuid;
using stringio::format;
using seqio::ChromosomeInstance;
using seqio::Locus;
using seqio::Nuc;
using seqio::TCoord;

namespace vario {

/* Mutation *
 *----------*/

Mutation::Mutation ()
: id(0), 
  is_snv(false),
  is_cnv(false)
{}

/* CopyNumberVariant *
 *-------------------*/

CopyNumberVariant::CopyNumberVariant ()
: id(0),
  is_wgd(false),
  is_deletion(false),
  is_chr_wide(false),
  is_telomeric(false),
  is_forward(true),
  len_rel(0.0),
  start_rel(0.0),
  ref_pos_begin(0),
  ref_pos_end(0),
  ref_chr("")
{}

ostream& operator<<(ostream& lhs, const CopyNumberVariant& cnv) {
  lhs << cnv.id;
  lhs << "\t" <<(cnv.is_wgd ? "WGD" : (cnv.is_deletion ? "DEL" : "CPY"));
  lhs << "\t" << (cnv.is_chr_wide ? "chr" : (cnv.is_telomeric ? "tel" : "foc"));
  lhs << "\t" << cnv.ref_chr;
  lhs << "\t" << cnv.start_rel;
  lhs << "\t" << cnv.is_forward ? "+" : "-";
  lhs << "\t" << cnv.len_rel;
  lhs << "\n";
  return lhs;
}
/* VariantAlleleCount *
 *--------------------*/

VariantAlleleCount::VariantAlleleCount ()
: num_tot(0), 
  num_alt(0)
{}

/*------------------------------------*/
/*           Utility methods          */
/*------------------------------------*/

unsigned assignSomaticMutationType(
  std::vector<Mutation>& vec_mutations,
  const double ratio_cnv,
  RandomNumberGenerator& rng)
{
  function<double()> random_double = rng.getRandomFunctionReal(0.0, 1.0);
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
    // var.reg_copy = 0; // TODO: can we do better than this?
    // var.chr_copy = 0; // TODO: can we do better than this?
    // TODO: at this point only SNVs are supported
    if (!var.isSnv()) {
      stringio::safeGetline(input, line);
      continue;
    }

    // for single-sample VCFs, set zygosity state
    if (num_samples == 1) {
      var.is_het = ( &var_cols[9][0] != &var_cols[9][2] );
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

void 
writeVcfHeader (
  ostream& out,
  vector<string> vec_seq_id,
  vector<TCoord> vec_seq_len,
  vector<string> vec_lbl
)
{
  // sanity check: do seq ids match seq lengths?
  assert ( vec_seq_id.size() == vec_seq_len.size() );

  out << "##fileformat=VCFv4.1" << endl;
  time_t timer = time(NULL);
  tm* t = localtime(&timer);
  out << format("##fileDate=%d-%d-%d", (1900+t->tm_year), t->tm_mon, t->tm_mday) << endl;
  out << "##source=CloniPhy v1.1.0" << endl;
  //out << "##reference=" << std::endl; # TODO: include ref filename
  vector<string> vec_ref_ids;
  for (size_t i=0; i<vec_seq_id.size(); i++) {
    string seq_id = vec_seq_id[i];
    TCoord seq_len = vec_seq_len[i];
    out << format("##contig=<ID=%s, length=%lu>", seq_id.c_str(), seq_len) << endl;
  }
  //out << "##phasing=complete" << endl;
  out << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" << endl;
  out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
  out << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << endl;
  
  // add column header
  out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for (auto lbl : vec_lbl) out << '\t' << lbl;
  out << endl;
}

/** Write multi-sample VCF (requires mutation matrix) */
void
writeVcfRecords (
  ostream& out,
  const vector<Variant>& vec_vars,
  const vector<vector<bool>>& mtx_mut, // binary mutation matrix (sample x var)
  const vector<int> vec_id_samples
)
{
  int num_samples = vec_id_samples.size();

  // defaults (introduce more params if necessary)
  string var_qual = "."; // QUAL column
  string var_info = format("NS=%d", num_samples); // INFO column
  string var_fmt  = "GT:GQ"; // FORMAT column
  short  gt_qual  = 60; // genotype quality (global)

  vector<Variant> sorted_vars = Variant::sortByPositionLex(vec_vars);
  for (Variant var : sorted_vars) {
    string ref = var.alleles[0];
    string alt = var.alleles[1];
    out << format("%s\t%d\t%s\t%s\t%s\t%s\tPASS\t%s\t%s",
                  var.chr.c_str(), (var.pos+1), var.id.c_str(),
                  ref.c_str(), alt.c_str(), var_qual.c_str(), 
                  var_info.c_str(), var_fmt.c_str());
    for (auto sid : vec_id_samples) {
      string genotype = mtx_mut[sid][var.idx_mutation] ? "0/1" : "0/0";
      out << format("\t%s:%d", genotype.c_str(), gt_qual);
    }
    out << endl;
  }
}

/** Write single-sample VCF (e.g., germline vars) */
void
writeVcfRecords (
  ostream& out,
  const vector<Variant>& vec_vars
)
{
  // defaults (introduce more params if necessary)
  string var_qual = "."; // QUAL column
  string var_info = "NS=1"; // INFO column
  string var_fmt  = "GT:GQ"; // FORMAT column
  short  gt_qual  = 60; // genotype quality (global)

  vector<Variant> sorted_vars = Variant::sortByPositionLex(vec_vars);
  for (Variant var : sorted_vars) {
    string ref = var.alleles[0];
    string alt = var.alleles[1];
    string var_geno = (var.is_het ? "0/1" : "1/1");
    out << format("%s\t%d\t%s\t%s\t%s\t%s\tPASS\t%s\t%s\t%s:%d",
                  var.chr.c_str(), (var.pos+1), var.id.c_str(),
                  ref.c_str(), alt.c_str(), var_qual.c_str(), var_info.c_str(), 
                  var_fmt.c_str(), var_geno.c_str(), gt_qual);
    out << endl;
  }
}


void writeVcf (
  const vector<shared_ptr<SeqRecord>>& seqs,
  const vector<Variant>& vars,
  const vector<int>& id_samples,
  const vector<string>& labels,
  const vector<vector<bool> >& mutMatrix,
  ostream& out)
{
  unsigned num_samples = id_samples.size();
  string var_qual = "."; // QUAL column
  string var_info = format("NS=%d", num_samples); // INFO column
  string var_fmt  = "GT:GQ"; // FORMAT column
  short  gt_qual  = 60; // genotype quality (global)

  // prepare ref seq ids and lengths
  vector<string> vec_ref_ids;
  vector<TCoord> vec_ref_len;
  for (auto rec : seqs) {
    vec_ref_ids.push_back(rec->id_ref);
    vec_ref_len.push_back(rec->seq.size());
  }

  // write header
  writeVcfHeader(out, vec_ref_ids, vec_ref_len, labels);

  // write variants
  vector<Variant> sorted_vars = Variant::sortByPositionLex(vars);
  for (Variant var : sorted_vars) {
    string ref = var.alleles[0];
    string alt = var.alleles[1];
    out << format("%s\t%d\t%s\t%s\t%s\t%s\tPASS\t%s\t%s",
                  var.chr.c_str(), (var.pos+1), var.id.c_str(),
                  ref.c_str(), alt.c_str(), var_qual.c_str(), 
                  var_info.c_str(), var_fmt.c_str());
    for (auto sid : id_samples) {
      string genotype = "";
      if (mutMatrix[sid][var.idx_mutation] == true)
        //genotype = (var.chr_copy==0 ? "1|0" : "0|1");
        genotype = (var.is_het ? "0/1" : "1/1");
      else
        //genotype = "0|0";
        genotype = "0/0";
      out << format("\t%s:%d", genotype.c_str(), gt_qual);
    }
    out << endl;
  }
}

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

unsigned 
writeVcf (
  const std::string filename,
  const GenomeReference& genome,
  const std::vector<Variant>& vars,
  const std::string label
)
{
  // defaults
  int num_samples = 1;
  string var_qual = "."; // QUAL column
  string var_info = format("NS=%d", num_samples); // INFO column
  string var_fmt  = "GT:GQ"; // FORMAT column
  short  gt_qual  = 60; // genotype quality (global)

  // open output file
  ofstream f_out;
  f_out.open(filename);
  
  // write header
  vector<string> vec_lbl = {label};
  writeVcfHeader(f_out, genome.vec_chr_id, genome.vec_chr_len, vec_lbl);

  // write variants
  writeVcfRecords(f_out, vars);

  // close output file
  f_out.close();  
}

long getRandomMutPos(
  function<int()> r_bucket,
  map<string, vector<long>> map_bucket_pos,
  short offset)
{
  return 0;
}

/** DEPRECATED! Functionality shifted to vario::VariantStore::generateGermlineVariants
    Generate variant loci in a given genome based on evolutionary model.
    Nucleotide substitution probabilities guide selection of loci. */
// vector<Variant> generateGermlineVariants (
//   const int num_variants,
//   const GenomeReference& genome,
//   GermlineSubstitutionModel& model,
//   RandomNumberGenerator<>& rng,
//   const bool inf_sites)
// {
//   vector<Variant> variants = vector<Variant>(num_variants);
//   boost::container::flat_set<int> var_pos; // keep track of variant positions
//   function<double()> random_float = rng.getRandomFunctionDouble(0.0, 1.0);
//   function<short()> random_copy = rng.getRandomFunctionInt(short(0), short(genome.ploidy-1));
//   random_selector<> selector(rng.generator); // used to pick random vector indices

//   // determine base mutation probs from model (marginal sums)
//   vector<double> p_i(4, 0);
//   for (int i=0; i<4; ++i) {
//     for (int j=0; j<4; ++j) {
//       p_i[i] += model.Qij[i][j];
//     }
//   }
//   function<int()> random_nuc_idx = rng.getRandomIndexWeighted(p_i);

//   unsigned long genome_len = genome.length; // haploid genome length
//   for (int i=0; i<num_variants; ++i) {
//     // pick random nucleotide bucket
//     int idx_bucket = random_nuc_idx();
//     // pick random position
//     long nuc_pos = selector(genome.nuc_pos[idx_bucket]);
//     if (inf_sites) {
//       while (binary_search(var_pos.begin(), var_pos.end(), nuc_pos)) {
// // TODO: check verbosity setting
// fprintf(stderr, "[INFO] Infinite sites assumption: locus %ld has been mutated before, picking another one...\n", nuc_pos);
//         nuc_pos = selector(genome.nuc_pos[idx_bucket]);
//       }
//       var_pos.insert(nuc_pos);
//     }
//     Locus loc = genome.getLocusByGlobalPos(nuc_pos);
//     // pick new nucleotide
//     short nuc_alt = evolution::MutateSite(idx_bucket, random_float, model);
//     Variant var;
//     var.id = str(format("g%d") % i);
//     var.chr = loc.id_ref;
//     //var.chr_copy = random_copy(); // TODO: deprecated!
//     var.rel_pos = double(nuc_pos-(var.chr_copy*genome_len))/genome_len;
//     var.reg_copy = 0; // TODO: deprecated!
//     var.pos = loc.start;
//     var.alleles.push_back(string(1, seqio::idx2nuc(idx_bucket)));
//     var.alleles.push_back(string(1, seqio::idx2nuc(nuc_alt)));
//     var.idx_mutation = i;
//     variants[i] = var;
//   }

//   return variants;
// }

/* !DEPRECATED! (use vario::VariantStore::generateGermlineVariants) */
/** Generate variant loci in a given genome based on evolutionary model.
    Nucleotide substitution probabilities guide selection of loci. 
vector<Variant> generateVariantsRandomPos(
  const int num_variants,
  const GenomeReference& genome,
  GermlineSubstitutionModel& model,
  RandomNumberGenerator& rng,
  const bool inf_sites)
{
  vector<Variant> variants = vector<Variant>(num_variants);
  boost::container::flat_set<int> var_pos; // keep track of variant positions
  function<double()> random_float = rng.getRandomFunctionReal(0.0, 1.0);
  function<long()> random_pos = rng.getRandomFunctionInt<long>(0, genome.length);
  random_selector<> selector(rng.generator); // used to pick random vector indices

  for (int i=0; i<num_variants; ++i) {
    // pick random position
    long nuc_pos = random_pos();
    if (inf_sites) {
      while (binary_search(var_pos.begin(), var_pos.end(), nuc_pos)) {
fprintf(stderr, "locus %ld has been mutated before, picking another one...\n", nuc_pos);
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
    //var.reg_copy = 0; // TODO: when implementing CNVs, use this property to indicate affected copy
    var.pos = nuc_pos;
    var.rel_pos = double(nuc_pos)/genome.length;
    var.alleles.push_back(string(1, seqio::idx2nuc(ref_nuc)));
    var.alleles.push_back(string(1, seqio::idx2nuc(nuc_alt)));
    variants[i] = var;
  }

  return variants;
} */

/* DEPRECATED! */
// void applyVariants(
//   GenomeReference &genome,
//   const vector<Variant> &variants)
// {
//   unsigned num_sequences = genome.num_records;
//   // generate lookup table for sequences
//   map<string,vector<unsigned>> chr2seq;
//   for (unsigned i=0; i<genome.records.size(); ++i) {
//     string id_ref = genome.records[i]->id_ref;
//     if (chr2seq.find(id_ref) == chr2seq.end()) {
//       chr2seq[id_ref] = vector<unsigned>();
//     }
//     chr2seq[id_ref].push_back(i);
//   }
// fprintf(stderr, "applying %lu variants...\n", variants.size());
//   // modify genome according to variant genotypes
//   for (Variant var : variants) {
//     auto it_chr_idx = chr2seq.find(var.chr);
//     // make sure ref seq exists in genome
//     if (it_chr_idx == chr2seq.end()) {
//       fprintf(stderr, "[WARN] sequence '%s' not contained in reference genome\n", var.chr.c_str());
//     }
//     else {
//       vector<unsigned> chr_idx = it_chr_idx->second;
//       // make sure mutated chr copy exists in genome
//       if (chr_idx.size() < var.chr_copy+1) {
//         fprintf(stderr, "[WARN] genome does not contain %d copies of CHR '%s'\n", var.chr_copy+1, var.chr.c_str());
//       }
//       else {
//         // apply variant to sequence
//         unsigned cidx = chr_idx[var.chr_copy];
//         genome.records[cidx]->seq[var.pos] = var.alleles[1][0]; // TODO: at the moment only SNVs are supported ("[0]" extracts the first character from the allel)
//       }
//     }
//   }
// }


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
