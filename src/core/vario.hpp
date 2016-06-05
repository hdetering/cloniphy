#ifndef VARIO_H
#define VARIO_H

#include "random.hpp"
#include "seqio.hpp"
#include "evolution.hpp"
#include <boost/function.hpp>
#include <fstream>
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
  short chr_copy;     /** affected copy of chromosome (0: mat, 1: pat, ...) */
  unsigned long pos;  /** reference basepair position */
  std::vector<std::string> alleles; /** observed alleles */
  unsigned idx_mutation; /** reference to mutation that gave rise to this variant */
  double rel_pos;     /** relative position in genome (use for sorting) */

  Variant();
  Variant(std::string id, std::string chr, unsigned long pos);
  ~Variant();

  bool operator< (const Variant&) const; /** make variants sortable */
  /** Sort variants by absolute position in genome. */
  static std::vector<Variant> sortByPosition(const std::vector<Variant>&);
  /** Returns true if this variant is a SNV, false otherwise. */
  bool isSnv();
};

/** VariantSets store a set of variants and summary statistics about them. */
struct VariantSet
{
  unsigned long num_variants = 0;
  std::vector<Variant> vec_variants; /** variants that belong to the set */
  /** summary statistics */
  double mat_freqs[4][4] = {
    { 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0, 0.0 }
  }; /** nucleotide substitution frequencies */

  VariantSet();
  ~VariantSet();
  long calculateSumstats();
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
  /** Apply single mutation to a given genomic sequence, */
  /** returning a Variant and Genotype object */
  void apply(
    Genome& genome,
    SubstitutionModel model,
    boost::function<double()>& random,
    Variant &var,
    Genotype &gt);
};

/** Generate random mutations out of thin air. */
std::vector<Mutation> generateMutations(
  const int num_mutations,
  boost::function<double()>&
);
/** Apply set of mutation to a given genomic sequence,
    returning a set of Variant objects
    (reference sequence is not changed) */
void applyMutations(
  const std::vector<Mutation> &,
  const Genome& genome,
  SubstitutionModel model,
  boost::function<double()>& random,
  std::vector<Variant> &variants);

/** Read VCF file and return list of variants. */
void readVcf(
  std::string vcf_filename,
  VariantSet& variants,
  std::vector<std::vector<Genotype> >& gtMatrix);
/** Read input stream with variants in VCF format and return list of variants. */
void readVcf(
  std::istream& vcf_filename,
  VariantSet& variants,
  std::vector<std::vector<Genotype> >& gtMatrix);
/** Generate VCF output from a reference sequence and a set of mutations.
    (multiple samples) */
void writeVcf(
  const std::vector<SeqRecord>& seqs,
  const std::vector<Variant>& vars,
  const std::vector<std::string>& labels,
  const std::vector<std::vector<short> >& mutMatrix,
  std::ostream&);
/** Generate VCF output from a reference sequence and a set of mutations.
    (single sample) */
void writeVcf(
  const std::vector<SeqRecord>& seqs,
  const std::vector<Variant>& vars,
  const std::string label,
  std::ostream&);

/** Generate variant loci in a given genome based on evolutionary model.
    Nucleotide substitution probabilities guide selection of loci. */
std::vector<Variant> generateVariants(
  const int num_variants,
  const Genome& genome,
  SubstitutionModel& model,
  RandomNumberGenerator<>&,
  const bool infinite_sites = false
);
/** Generate variant loci in a given genome based on evolutionary model.
    Loci are selected randomly. */
std::vector<Variant> generateVariantsRandomPos(
  const int num_variants,
  const Genome& genome,
  SubstitutionModel& model,
  RandomNumberGenerator<>&,
  const bool infinite_sites = false
);
/** Apply variants to a given reference sequence */
void applyVariants(
  Genome&,
  const std::vector<Variant>&,
  const std::vector<Genotype>&);
/** Apply variants to a given reference sequence. streaming modified genome to output. */
void applyVariantsStream(
  const Genome &ref_genome,
  const std::vector<Mutation> &mutations,
  const std::vector<Variant> &variants,
  std::ostream &outstream,
  short len_line = 60);
} /* namespace vario */

#endif /* VARIO_H */
