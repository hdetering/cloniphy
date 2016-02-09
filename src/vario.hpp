#ifndef VARIO_H
#define VARIO_H

#include "seqio.hpp"
#include "evolution.hpp"
#include <boost/function.hpp>
#include <iostream>
#include <string>
#include <vector>

using seqio::Genome;
using seqio::SeqRecord;
using evolution::SubstitutionModel;

namespace vario {

/** Variants represent variable sites in nucleotide sequences */
struct Variant
{
  std::string id;     /** unique identifier */
  std::string chr;    /** reference chromosome id */
  unsigned long pos;  /** reference basepair position */
  std::vector<std::string> alleles; /** observed alleles */
  unsigned idx_mutation; /** reference to mutation that gave rise to this variant */
  double rel_pos;     /** relative position in genome (use for sorting) */

  Variant();
  ~Variant();

  bool operator< (const Variant&) const; /** make variants sortable */
  /** Sort variants by absolute position in genome. */
  static std::vector<Variant> sortByPosition(const std::vector<Variant>&);
  /** Returns true if this variant is a SNV, false otherwise. */
  bool isSnv();
};

struct Genotype
{
  std::string id_variant; /** the variant this genotype refers to */
  short maternal; /** allele on maternal strand */
  short paternal; /** allele on paternal strand */
};

/** Mutations specifiy a modification of a sequence.
 *  TODO: Do I really belong here?
 */
struct Mutation
{
  unsigned id;   /** unique identifier */
  double relPos; /** relative position in genome [0..1] */
  short copy;    /** which chromosome copy (0:maternal, 1:paternal) */

  Mutation();
  //Mutation(char ref, char alt); /** c'tor converts alleles to nucleotide shift. */

  bool operator< (const Mutation&) const; /** make mutations sortable */
  static std::vector<Mutation> sortByPosition(const std::vector<Mutation>&);
  /** Apply mutation to a given genomic sequence, */
  /** returning a Variant and Genotype object */
  void apply(
    Genome& genome,
    SubstitutionModel model,
    boost::function<float()>& random,
    Variant &var,
    Genotype &gt);
};

/** Generate random mutations out of thin air. */
std::vector<Mutation> generateMutations(
  const int num_mutations,
  boost::function<float()>&
);
/** Read VCF file and return list of variants. */
void readVcf(
  std::string vcf_filename,
  std::vector<Variant>& variants,
  std::vector<std::vector<Genotype> >& gtMatrix);
/** Read input stream with variants in VCF format and return list of variants. */
void readVcf(
  std::istream& vcf_filename,
  std::vector<Variant>& variants,
  std::vector<std::vector<Genotype> >& gtMatrix);
/** Generate VCF output from a reference sequence and a set of mutations.  */
void writeVcf(
  const std::vector<SeqRecord>&,
  const std::vector<Variant>&,
  const std::vector<std::string>&,
  const std::vector<std::vector<short> >&,
  std::ostream&);
/** Apply variants to a given reference sequence */
void applyVariants(
  Genome&,
  const std::vector<Variant>&,
  const std::vector<Genotype>&);

} /* namespace vario */

#endif /* VARIO_H */
