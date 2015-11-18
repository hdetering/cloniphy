#include "vario.hpp"
#include <algorithm>
#include <boost/format.hpp>
#include <ctime>
#include <string>

using boost::format;
using std::vector;
using std::string;

bool Mutation::operator< (const Mutation &other) const {
  return absPos < other.absPos;
}

vector<Mutation> Mutation::sortByPosition(const vector<Mutation> &mutations) {
  std::vector<Mutation> mutationsCopy = mutations;
  std::sort(mutationsCopy.begin(), mutationsCopy.end());
  return mutationsCopy;
}

void VarIO::writeVcf(const vector<SeqRecord> &seqs, const vector<Mutation> &muts, const vector<vector<short> > &mutMatrix, std::ostream &out) {
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
  unsigned long cum_len = rec.seq.length();
  for (vector<Mutation>::iterator mut=sorted_muts.begin(); mut!=sorted_muts.end(); ++mut) {
    unsigned long p = mut->absPos;
    // check if we passed the end of the current sequence
    if (p >= cum_len) {
      rec = seqs[++idx_seq];
      cum_len += rec.seq.length();
    }
    char ref = rec.seq[p];
    char alt = SeqIO::shiftNucleotide(ref, mut->offset);
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
