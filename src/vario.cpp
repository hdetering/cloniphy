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

namespace vario {

Mutation::Mutation() {
  this->id = 0;
  this->absPos = 0;
  this->offset = 0;
  this->copy = 0;
}

Mutation::Mutation(char ref, char alt) {
  short ref_nuc = seqio::charToNuc(ref);
  short alt_nuc = seqio::charToNuc(alt);
  this->offset =  ((alt_nuc-ref_nuc) % 4); // TODO: this could be more generic (get rid of the hard-coded 4)
}

bool Mutation::operator< (const Mutation &other) const {
  return absPos < other.absPos;
}

vector<Mutation> Mutation::sortByPosition(const vector<Mutation> &mutations) {
  vector<Mutation> mutationsCopy = mutations;
  sort(mutationsCopy.begin(), mutationsCopy.end());
  return mutationsCopy;
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

vector<Mutation> generateMutations(const int num_mutations, unsigned long ref_len, boost::function<float()>& random) {

  vector<Mutation> mutations(num_mutations);
  std::set<unsigned long> mutPositions; // remember mutated positions (enforce infinite sites model)
  unsigned i=0;

  for (vector<Mutation>::iterator m=mutations.begin(); m!=mutations.end(); ++m) {
    float rel_pos = random();
    unsigned long abs_pos = rel_pos * ref_len;
    // enforce infinite sites (no position can be mutated more than once)
    while (mutPositions.count(abs_pos)!=0) {
      abs_pos = (abs_pos+1) % ref_len; // wrap around at end of reference
    }
    short offset = (random()*3)+1;
    short copy = random()*2;
//fprintf(stderr, "<Mutation(id=%u;absPos=%ld,offset=%d,copy=%d)>\n", i, abs_pos, offset, copy);
    m->id = i++;
    m->absPos = abs_pos;
    m->offset = offset;
    m->copy = copy;
  }

  return mutations;
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
      Genotype gt = {a_allele, b_allele};
      //Genotype *gt = new Genotype();
      //gt->maternal = a_allele;
      //gt->paternal = b_allele;
      gtMatrix[i].push_back(gt);
    }
    stringio::safeGetline(input, line);
  }
}

void writeVcf(const vector<SeqRecord> &seqs, const vector<Mutation> &muts, const vector<vector<short> > &mutMatrix, std::ostream &out) {
  unsigned num_samples = mutMatrix.size();
  short  var_qual = 40; // QUAL column
  string var_info = (format("NS=%d") % num_samples).str(); // INFO column
  string var_fmt  = "GT:GQ"; // FORMAT column
  short  gt_qual  = 60; // genotype quality (global)

  // write header
  out << "##fileformat=VCFv4.1" << std::endl;
  time_t timer = time(NULL);
  tm* t = localtime(&timer);
  out << "##fileDate=" << (1900+t->tm_year) << t->tm_mon << t->tm_mday << std::endl;
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
    out << format("\tclone%02d") % i;
  }
  out << std::endl;

  // write variants
  vector<Mutation> sorted_muts = Mutation::sortByPosition(muts);
  unsigned idx_seq = 0;
  SeqRecord rec = seqs[0];
  unsigned long seq_len = rec.seq.length();
  unsigned long cum_len = 0;
  for (vector<Mutation>::iterator mut=sorted_muts.begin(); mut!=sorted_muts.end(); ++mut) {
    unsigned long p = mut->absPos;
    // check if we passed the end of the current sequence
    while (p >= (cum_len+seq_len)) {
      rec = seqs[++idx_seq];
      seq_len = rec.seq.length();
      cum_len += seq_len;
    }
    char ref = rec.seq[p-cum_len];
    char alt = seqio::shiftNucleotide(ref, mut->offset);
    out << format("%s\t%d\t%d\t%s\t%s\t%d\tPASS\t%d\t%s") % rec.id % p % mut->id % ref % alt % var_qual % var_info % var_fmt;
    for (unsigned i=0; i<num_samples; ++i) {
      string genotype = "";
      if (mutMatrix[i][mut->id] != 0) {
        genotype = (mut->copy==0 ? "1|0" : "0|1"); }
      else {
        genotype = "0|0"; }
      out << format("\t%s:%d") % genotype % gt_qual;
    }
    out << std::endl;
  }
}

void applyVariants(vector<SeqRecord>& sequences, const vector<Variant>& variants, const vector<Genotype>& genotypes) {
  unsigned num_sequences = sequences.size();
  // generate lookup table for sequences
  map<string,unsigned> chr2seq;
  for (unsigned i=0; i<num_sequences; ++i) {
    chr2seq[sequences[i].id] = i;
  }
  // duplicate reference
  for (unsigned i=0; i<num_sequences; ++i) {
    SeqRecord orig = sequences[i];
    SeqRecord *dupl = new SeqRecord(orig.id, orig.description, orig.seq);
    sequences.push_back(*dupl);
    sequences[i].id += "_m";
    sequences[num_sequences+i].id += "_p";
  }
fprintf(stderr, "variants: %lu\n", variants.size());
fprintf(stderr, "genotypes: %lu\n", genotypes.size());
  // modify reference according to variant genotypes
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
        sequences[chr_idx].seq[var.pos-1] = var.alleles[gt.maternal][0]; // TODO: at the moment only SNVs are supported ("[0]" extracts the first character from the allel)
      }
      if (gt.paternal>0) {
//fprintf(stderr, "\t%lu paternal: '%s' -> '%s'\n", var.pos, var.alleles[0].c_str(), var.alleles[gt.paternal].c_str());
        sequences[num_sequences+chr_idx].seq[var.pos-1] = var.alleles[gt.paternal][0];
      }
    }
    else {
      fprintf(stderr, "[WARN] sequence '%s' not contained in reference genome\n", var.chr.c_str());
    }
  }
}

} /* namespace vario */
