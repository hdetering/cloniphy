#ifndef VARIO_H
#define VARIO_H

#include "random.hpp"
#include "seqio.hpp"
#include "stringio.hpp"
#include "evolution.hpp"
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using seqio::Genome;
using seqio::SeqRecord;
using evolution::GermlineSubstitutionModel;
using evolution::SomaticSubstitutionModel;

/**
 * Classes and methods for input/output and simulation of variants.
 */
namespace vario {

/** Variants represent variable sites in nucleotide sequences */
struct Variant
{
  std::string id;     /** unique identifier */
  std::string chr;    /** reference chromosome id */
  short chr_copy;     /** chromosome phase (0: maternal, 1: paternal, ...[hyperdiploid]) */
  short reg_copy;     /** affected copy of chr region (0: original, 1: first copy, ...) */
  unsigned long pos;  /** reference basepair position */
  std::vector<std::string> alleles; /** observed alleles */
  unsigned idx_mutation; /** reference to mutation that gave rise to this variant */
  double rel_pos;     /** relative position in genome (use for sorting) */
  bool is_somatic;    /** is this Variant somatic or germline? (different output channels) */

  Variant();
  Variant(std::string id, std::string chr, unsigned long pos);
  ~Variant();

  bool operator< (const Variant&) const; /** make variants sortable */
  /** Sort variants by absolute position in genome. */
  static std::vector<Variant> sortByPosition(const std::vector<Variant>&);
  /** Sort variants, first lexicographically by CHR, then by position in CHR. */
  static std::vector<Variant> sortByPositionLex(const std::vector<Variant>&);
  /** Sort variants by absolute position in genome, taking chomosome copies into account. */
  static std::vector<Variant> sortByPositionPoly(const std::vector<Variant>&);
  /** Sort variants by position in reference genome. */
  static std::vector<Variant> sortByPositionRef(const std::vector<Variant>&);
  /** Returns true if this variant is a SNV, false otherwise. */
  bool isSnv();
};

/** VariantSets store a set of variants and summary statistics about them. */
struct VariantSet
{
  unsigned long num_variants = 0;
  std::vector<Variant> vec_variants; /** all variants that belong to the set */
  // TODO: a position must be able to store more than one variable! (-> CNVs, homoplasy)
  std::map<std::string, std::map<unsigned long, std::vector<Variant>>> map_chr2pos2var; /** variants stored by chromosome id */

  /** summary statistics */
  double mat_freqs[4][4] = {
    { 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.0, 0.0 }
  }; /** nucleotide substitution frequencies */

  VariantSet();
  VariantSet(std::vector<Variant> variants);
  ~VariantSet();

  /** Index variants by chromosome and position. */
  long indexVariants();
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
    GermlineSubstitutionModel model,
    std::function<double()>& random,
    Variant &var,
    Genotype &gt);
};

/** Generate random mutations out of thin air. */
std::vector<Mutation> generateMutations(
  const int num_mutations,
  std::function<double()>&
);
/** Apply set of mutations to a given genomic sequence,
    returning a set of Variant objects
    (reference sequence is not changed) */
void applyMutations(
  const std::vector<Mutation> &,
  const Genome& genome,
  GermlineSubstitutionModel model,
  std::function<double()>& random,
  std::vector<Variant> &variants);

/** Read VCF file and return list of variants. */
void readVcf(
  std::string fn_vcf,
  VariantSet& variants,
  std::map<std::string, std::vector<Genotype> >& gtMatrix);
/** Read input stream with variants in VCF format and return list of variants. */
void readVcf(
  std::istream& fs_vcf,
  VariantSet& variants,
  std::map<std::string, std::vector<Genotype> >& gtMatrix);
/** Generate VCF output for a reference genome and a set of mutations.
    (multiple samples) */
void writeVcf(
  const std::vector<SeqRecord>& seqs,
  const std::vector<Variant>& vars,
  const std::vector<int>& id_samples,
  const std::vector<std::string>& labels,
  const std::vector<std::vector<bool> >& mutMatrix,
  std::ostream&);
/** Generate VCF output from a reference genome and a set of variants.
    (single sample) */
void writeVcf(
  const std::vector<SeqRecord>& seqs,
  const std::vector<Variant>& vars,
  const std::string label,
  const std::string filename);

/** Read mutation map (clone x mutation) from a CSV file. */
int readMutMapFromCSV(
  std::map<std::string, std::vector<bool>> &mm,
  const std::string &filename
);

/** Generate variant loci in a given genome based on evolutionary model.
    Nucleotide substitution probabilities guide selection of loci. */
std::vector<Variant> generateGermlineVariants(
  const int num_variants,
  const Genome& genome,
  GermlineSubstitutionModel& model,
  RandomNumberGenerator<>&,
  const bool infinite_sites = false
);
/** Generate variant loci in a given genome based on somatic mutation model.
    Use context-dependent mutation signature to select loci. */
std::vector<Variant> generateSomaticVariants(
  const int num_variants,
  const Genome& genome,
  SomaticSubstitutionModel& model,
  RandomNumberGenerator<>&,
  const bool infinite_sites = false
);
/** Generate variant loci in a given genome based on evolutionary model.
    Loci are selected randomly. */
std::vector<Variant> generateVariantsRandomPos(
  const int num_variants,
  const Genome& genome,
  GermlineSubstitutionModel& model,
  RandomNumberGenerator<>&,
  const bool infinite_sites = false
);
/** Apply variants to a given reference sequence */
void applyVariants(
  Genome&,
  const std::vector<Variant>&);
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
