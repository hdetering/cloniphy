#ifndef VARIO_H
#define VARIO_H

#include "seqio.hpp"
#include <boost/function.hpp>
#include <iostream>
#include <string>
#include <vector>

using seqio::SeqRecord;

namespace vario {

/** Mutations specifiy a modification of a sequence.
 *  TODO: Do I really belong here?
 */
struct Mutation
{
  unsigned id;          /** unique identifier */
  unsigned long absPos; /** absolute bp position */
  short offset;         /** shifts the ancestral genotype to a new one */
  short copy;           /** which chromosome copy (0:maternal, 1:paternal) */

  Mutation();
  Mutation(char ref, char alt); /** c'tor converts alleles to nucleotide shift. */

  bool operator< (const Mutation&) const; /** make mutations sortable */
  static std::vector<Mutation> sortByPosition(const std::vector<Mutation>&);
};

/** Variants represent variable sites in nucleotide sequences */
struct Variant
{
  std::string id;    /** unique identifier */
  std::string chr;   /** reference chromosome id */
  unsigned long pos; /** reference basepair position */
  std::vector<std::string> alleles; /** observed alleles */

  /** Returns true if this variant is a SNV, false otherwise. */
  bool isSnv();
};

struct Genotype
{
  short maternal; /** allele on maternal strand */
  short paternal; /** allele on paternal strand */
};

/** Generate random mutations out of thin air. */
std::vector<Mutation> generateMutations(const int num_mutations, unsigned long ref_len, boost::function<float()>& random);
/** Read VCF file and return list of variants. */
void readVcf(std::string vcf_filename, std::vector<Variant>& variants, std::vector<std::vector<Genotype> >& gtMatrix);
/** Read input stream with variants in VCF format and return list of variants. */
void readVcf(std::istream& vcf_filename, std::vector<Variant>& variants, std::vector<std::vector<Genotype> >& gtMatrix);
/** Generate VCF output from a reference sequence and a set of mutations.  */
void writeVcf(const std::vector<SeqRecord>&, const std::vector<Mutation>&, const std::vector<std::vector<short> >&, std::ostream&);
/** Apply variants to a given reference sequence */
void applyVariants(std::vector<SeqRecord>&, const std::vector<Variant>&, const std::vector<Genotype>&);

} /* namespace vario */

#endif /* VARIO_H */
