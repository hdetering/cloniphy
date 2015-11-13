#include "vario.hpp"
#include <algorithm>
#include <ctime>

bool Mutation::operator< (const Mutation &other) const {
  return absPos < other.absPos;
}

std::vector<Mutation> Mutation::sortByPosition(const std::vector<Mutation> &mutations) {
  std::vector<Mutation> mutationsCopy = mutations;
  std::sort(mutationsCopy.begin(), mutationsCopy.end());
  return mutationsCopy;
}

void VarIO::writeVcf(const std::vector<SeqRecord> &seqs, const std::vector<Mutation> &muts, std::ostream &out) {
  // write header
  out << "##fileformat=VCFv4.1" << std::endl;
  time_t timer = time(NULL);
  tm* t = localtime(&timer);
  out << "##fileDate=" << (1900+t->tm_year) << t->tm_mon << t->tm_mday << std::endl;
  out << "##source=CloniPhy v0.01" << std::endl;
  //out << "##reference=" << std::endl; # TODO: include ref filename
  out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;

  // write variants
  std::vector<Mutation> sorted_muts = Mutation::sortByPosition(muts);
  unsigned idx_seq = 0;
  SeqRecord rec = seqs[0];
  unsigned long cum_len = rec.seq.length();
  for (std::vector<Mutation>::iterator mut=sorted_muts.begin(); mut!=sorted_muts.end(); ++mut) {
    unsigned long p = mut->absPos;
    // check if we passed the end of the current sequence
    if (p >= cum_len) {
      rec = seqs[++idx_seq];
      cum_len += rec.seq.length();
    }
    char ref = rec.seq[p];
    char alt = SeqIO::shiftNucleotide(ref, mut->offset);
    out << rec.id<<"\t"<<p<<"\t"<<"0"<<"\t"<<ref<<"\t"<<alt<<"\t"<<"40\tPASS\t"<<std::endl;
  }
}
