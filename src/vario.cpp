#include "vario.hpp"
#include <algorithm>
#include <boost/format.hpp>
#include <ctime>
#include <fstream>
#include <map>
#include <set>
#include <stdio.h>

using namespace std;
using boost::format;
using seqio::Locus;
using seqio::Nuc;

namespace vario {

Mutation::Mutation() {
  this->id = 0;
  this->relPos = 0.0;
  this->copy = 0;
}

/*
Mutation::Mutation(char ref, char alt) {
  short ref_nuc = seqio::charToNuc(ref);
  short alt_nuc = seqio::charToNuc(alt);
  this->offset =  ((alt_nuc-ref_nuc) % 4); // TODO: this could be more generic (get rid of the hard-coded 4)
}*/

bool Mutation::operator< (const Mutation &other) const {
  return (relPos + copy) < (other.relPos + other.copy);
}

vector<Mutation> Mutation::sortByPosition(const vector<Mutation> &mutations) {
  vector<Mutation> mutationsCopy = mutations;
  sort(mutationsCopy.begin(), mutationsCopy.end());
  return mutationsCopy;
}

void Mutation::apply(
       Genome &genome,
       SubstitutionModel model,
       boost::function<float()> &rng,
       Variant &var,
       Genotype &gt)
{
  Locus loc = genome.getAbsoluteLocusMasked(this->relPos);
  Nuc nuc_ref = seqio::charToNuc(genome.records[loc.idx_record].seq[loc.start]);
  Nuc nuc_alt = (Nuc)(model.MutateNucleotide((short)nuc_ref, rng));

  // initialize variant
  var.id = str(boost::format("m%u") % this->id);
  var.chr = genome.records[loc.idx_record].id_ref;
  var.pos = loc.start;
  var.alleles.push_back(seqio::nucToString(nuc_ref));
  var.alleles.push_back(seqio::nucToString(nuc_alt));
  var.idx_mutation = this->id;
  var.rel_pos = this->relPos;

  // derive genotype from mutation
  gt.id_variant = var.id;
  gt.maternal = (this->copy==0 ? 1 : 0);
  gt.paternal = (this->copy==1 ? 1 : 0);

  // modify genomic sequence
  unsigned targetSeqIndex = loc.idx_record + (genome.num_records * this->copy);
  genome.records[targetSeqIndex].seq[loc.start] = seqio::nucToChar(nuc_alt);

fprintf(stderr, "Mutated '%s:%u' (%s>%s)\n", genome.records[targetSeqIndex].id.c_str(), loc.start, var.alleles[0].c_str(), var.alleles[1].c_str());
}

Variant::Variant() : id(""), chr(""), pos(0), alleles(0), idx_mutation(0), rel_pos(0.0) {}
Variant::~Variant() {}

bool Variant::operator< (const Variant &other) const {
  return rel_pos < other.rel_pos;
}

vector<Variant> Variant::sortByPosition(const vector<Variant> &variants) {
  vector<Variant> variantsCopy = variants;
  sort(variantsCopy.begin(), variantsCopy.end());
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

vector<Mutation> generateMutations(const int num_mutations, boost::function<float()>& random) {
  vector<Mutation> mutations(num_mutations);
  unsigned i=0;

  for (vector<Mutation>::iterator m=mutations.begin(); m!=mutations.end(); ++m) {
    float rel_pos = random();
    short copy = random()*2;
//fprintf(stderr, "<Mutation(id=%u;relPos=%f,copy=%d)>\n", i, rel_pos, copy);
    m->id = i++;
    m->relPos = rel_pos;
    m->copy = copy;
  }

  return mutations;
}

/** Apply set of mutation to a given genomic sequence,
    returning a set of Variant objects
    (actual reference sequence is not changed) */
void applyMutations(
  const std::vector<Mutation> & mutations,
  const Genome& genome,
  SubstitutionModel model,
  boost::function<float()>& random,
  vector<Variant> &variants)
{
  for (vector<Mutation>::const_iterator m=mutations.begin(); m!=mutations.end(); ++m) {
    Locus loc = genome.getAbsoluteLocusMasked(m->relPos);
    Nuc nuc_ref = seqio::charToNuc(genome.records[loc.idx_record].seq[loc.start]);
    Nuc nuc_alt = (Nuc)(model.MutateNucleotide((short)nuc_ref, random));

    // initialize variant
    Variant var = *(new Variant());
    var.id = str(boost::format("m%u") % m->id);
    var.chr = genome.records[loc.idx_record].id_ref;
    var.pos = loc.start;
    var.alleles.push_back(seqio::nucToString(nuc_ref));
    var.alleles.push_back(seqio::nucToString(nuc_alt));
    var.idx_mutation = m->id;
    var.rel_pos = m->relPos;

    // append variant to output
    variants.push_back(var);
  }
}

void readVcf(string vcf_filename, vector<Variant>& variants, vector<vector<Genotype> > &gtMatrix) {
  ifstream f_vcf;
  f_vcf.open(vcf_filename.c_str(), ios::in);
  readVcf(f_vcf, variants, gtMatrix);
  f_vcf.close();
}

void readVcf(std::istream &input, vector<Variant>& variants, vector<vector<Genotype> > &gtMatrix) {
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
  num_samples = hdr_cols.size()-9; // samples start at column 10
  gtMatrix = vector<vector<Genotype> >(num_samples);

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
    var.pos = atol(var_cols[1].c_str());
    var.id  = var_cols[2];
    var.alleles.push_back(var_cols[3]);
    vector<string> alt = stringio::split(var_cols[4], ',');
    var.alleles.insert(var.alleles.end(), alt.begin(), alt.end());
    // TODO: at this point only SNVs are supported
    if (!var.isSnv()) {
      stringio::safeGetline(input, line);
      continue;
    }

    variants.push_back(var);
    var_idx++;
//fprintf(stderr, "read from VCF: <Variant(idx=%u,id=%s,pos=%lu,num_alleles=%lu)> ", var_idx, var->id.c_str(), var->pos, var->alt.size()+1);
    for (unsigned i=0; i<num_samples; ++i) {
      string genotype = var_cols[9+i];
      short a_allele = atoi(&genotype[0]);
      short b_allele = atoi(&genotype[2]); // TODO: this will fail for >9 alleles.
//fprintf(stderr, "Genotype for sample %u: %d, %d\n", i, a_allele, b_allele);
      Genotype gt = { var.id, a_allele, b_allele };
      gtMatrix[i].push_back(gt);
    }
    stringio::safeGetline(input, line);
  }
}

void writeVcf(const vector<SeqRecord> &seqs, const vector<Variant> &vars, const vector<string> &labels, const vector<vector<short> > &mutMatrix, std::ostream &out) {
  unsigned num_samples = mutMatrix.size();
  short  var_qual = 40; // QUAL column
  string var_info = (format("NS=%d") % num_samples).str(); // INFO column
  string var_fmt  = "GT:GQ"; // FORMAT column
  short  gt_qual  = 60; // genotype quality (global)

  // write header
  out << "##fileformat=VCFv4.1" << std::endl;
  time_t timer = time(NULL);
  tm* t = localtime(&timer);
  out << boost::format("##fileDate=%d-%d-%d") % (1900+t->tm_year) % t->tm_mon % t->tm_mday << std::endl;
  out << "##source=CloniPhy v0.01" << std::endl;
  //out << "##reference=" << std::endl; # TODO: include ref filename
  for (vector<SeqRecord>::const_iterator ref=seqs.begin(); ref!=seqs.end(); ++ref) {
    out << format("##contig=<ID=%d, length=%u>") % ref->id % ref->seq.size() << std::endl;
  }
  out << "##phasing=complete" << std::endl;
  out << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">" << std::endl;
  out << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
  out << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << std::endl;
  out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  out << "\thealthy";
  for (unsigned i=1; i<num_samples; ++i) {
    out << '\t' << labels[i];
    //out << format("\tclone%02d") % i;
  }
  out << std::endl;

  // write variants
  vector<Variant> sorted_vars = Variant::sortByPosition(vars);
  for (vector<Variant>::iterator var=sorted_vars.begin(); var!=sorted_vars.end(); ++var) {
    string ref = var->alleles[0];
    string alt = var->alleles[1];
    out << format("%s\t%d\t%d\t%s\t%s\t%d\tPASS\t%d\t%s")
                  % (var->chr) % (var->pos+1) % (var->idx_mutation)
                  % ref % alt % var_qual % var_info % var_fmt;
    for (unsigned i=0; i<num_samples; ++i) {
      string genotype = "";
      short copy = mutMatrix[i][var->idx_mutation];
      if (copy != 0) {
        genotype = (copy==1 ? "1|0" : "0|1"); }
      else {
        genotype = "0|0"; }
      out << format("\t%s:%d") % genotype % gt_qual;
    }
    out << std::endl;
  }
}

void applyVariants(Genome &genome, const vector<Variant> &variants, const vector<Genotype> &genotypes) {
  unsigned num_sequences = genome.num_records;
  // generate lookup table for sequences
  map<string,unsigned> chr2seq;
  for (unsigned i=0; i<num_sequences; ++i) {
    chr2seq[genome.records[i].id_ref] = i;
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
        genome.records[chr_idx].seq[var.pos-1] = var.alleles[gt.maternal][0]; // TODO: at the moment only SNVs are supported ("[0]" extracts the first character from the allel)
      }
      if (gt.paternal>0) {
//fprintf(stderr, "\t%lu paternal: '%s' -> '%s'\n", var.pos, var.alleles[0].c_str(), var.alleles[gt.paternal].c_str());
        genome.records[num_sequences+chr_idx].seq[var.pos-1] = var.alleles[gt.paternal][0];
      }
    }
    else {
      fprintf(stderr, "[WARN] sequence '%s' not contained in reference genome\n", var.chr.c_str());
      return;
    }
  }
}

void applyVariantsStream(
  const Genome &g,
  const vector<Mutation> &mutations,
  const vector<Variant> &variants,
  ostream &outstream,
  short len_line)
{
  // sort mutations according to (diploid) genomic position
  vector<Mutation> vec_mut_sorted = Mutation::sortByPosition(mutations);
  unsigned idx_mut = 0;
  unsigned idx_chr = 0;
  unsigned idx_nuc = 0;
  short idx_char = 0;
  bool chr_mutated = false;
  Variant next_var = variants[vec_mut_sorted[idx_mut].id];
  for (idx_chr=0; idx_chr<g.records.size(); ++idx_chr) {
    chr_mutated = (g.records[idx_chr].id_ref == next_var.chr); // is this chromosome mutated?
    chr_mutated = chr_mutated && (g.records[idx_chr].copy == vec_mut_sorted[idx_mut].copy); // is this sequence the right copy?
    // print header
    outstream << ">" << g.records[idx_chr].id << endl;
    idx_nuc = 0;
    for (string::const_iterator nuc=g.records[idx_chr].seq.begin();
         nuc!=g.records[idx_chr].seq.end(); ++nuc) {
      if (chr_mutated && idx_nuc++ == next_var.pos) {
        outstream << next_var.alleles[1];
        if (++idx_mut<vec_mut_sorted.size())
          next_var = variants[vec_mut_sorted[idx_mut].id];
        else // print reference nucleotide
          chr_mutated = false;
      }
      else
        outstream << *nuc;
      if (++idx_char == len_line) { // enforce fixed line width
        outstream << endl;
        idx_char = 0;
      }
    }
    if (idx_char != 0) // end of sequence
      outstream << endl;
  }
}

} /* namespace vario */
