#ifndef VARIO_H
#define VARIO_H

#include "random.hpp"
#include "seqio.hpp"
#include "seqio/ChromosomeInstance.hpp"
#include "seqio/GenomeReference.hpp"
#include "seqio/GenomeInstance.hpp"
#include "stringio.hpp"
#include "evolution.hpp"
#include "vario/Variant.hpp"
#include "vario/VariantSet.hpp"
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

using seqio::GenomeReference;
using seqio::GenomeInstance;
using seqio::SeqRecord;
using seqio::SegmentCopy;
using evolution::GermlineSubstitutionModel;
using evolution::SomaticSubstitutionModel;
using evolution::SomaticCnvModel;

/**
 * Classes and methods for input/output and simulation of variants.
 */
namespace vario {

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
  /** unique identifier */
  int id;
  /** true: this Mutation is a single-nucleotide variant */
  bool is_snv;
  /** true: this Mutation is a copy-number variant */
  bool is_cnv;

  Mutation();

  // bool operator< (const Mutation&) const; /** make mutations sortable */
  // static std::vector<Mutation> sortByPosition(const std::vector<Mutation>&);
};

/** CopyNumberVariants encapsulate CNV events with all of their properties */
struct CopyNumberVariant
{
  unsigned    id;              /** unique identifier */
  bool        is_wgd;          /** true: CNV is Whole Genome Duplication */
  bool        is_deletion;     /** true: CNV is deletion event */
  bool        is_chr_wide;     /** true: event affects whole chromosome */
  bool        is_telomeric;    /** true: event coordinates include chromosome end */
  bool        is_forward;      /** true: event at 3' side of start_rel; false: at 5' end */
  double      len_rel;         /** length of affected region (fraction of chromsome length) */
  double      start_rel;       /** start position of event (fraction of chromosome length) */
  unsigned    ref_pos_begin;   /** start coordinate (in reference chr) */
  unsigned    ref_pos_end;     /** end coordinate (in reference chr) */
  std::string ref_chr;         /** affected chromosome (reference ID) */

  /** default c'tor */
  CopyNumberVariant();
};

/** Encapsulates the allele counts (reference, alternative) for a given variant and genome. */
struct VariantAlleleCount {
  //int idx_var; // this should rather be an index under which to store VariantAlleleCount elements.
  /** Total number of copies. */
  short num_tot;
  /** Number of alternative alleles. */
  short num_alt;

  /** Default c'tor. */
  VariantAlleleCount ();
};

//typedef std::map<seqio::TCoord, std::tuple<vario::Variant, short, short>> TMapPosVaf;
typedef std::map<seqio::TCoord, VariantAlleleCount> TMapPosVaf;
typedef std::map<std::string, TMapPosVaf> TMapChrPosVaf;

/** Inititalize a list of mutations, assigning a type (single-nucleotide vs. copy-number).
 *  \param vec_mutations list of (uninitialized) mutation objects
 *  \param ratio_cnv fraction of mutations that should be assigned CNV type
 *  \param rng RandomNumberGenerator object (reproducability)
 *  \return number of CNV mutations
 */
unsigned assignSomaticMutationType(
  std::vector<Mutation>& vec_mutations,
  const double ratio_cnv,
  RandomNumberGenerator& rng);

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
  const std::vector<std::shared_ptr<SeqRecord>>& seqs,
  const std::vector<Variant>& vars,
  const std::vector<int>& id_samples,
  const std::vector<std::string>& labels,
  const std::vector<std::vector<bool> >& mutMatrix,
  std::ostream&);

/** Generate VCF output from a reference genome and a set of variants (single sample).
  * \param seqs      Reference sequences (IDs and lengths are reported in VCF header).
  * \param vars      Variants to be output.
  * \param label     Header for the genotype column.
  * \param filename  Output file name.
  */
void 
writeVcf (
  const std::vector<std::shared_ptr<SeqRecord>>& seqs,
  const std::vector<Variant>& vars,
  const std::string label,
  const std::string filename
);

/** Generate VCF output from a reference genome and a set of variants (single sample).
  * \param filename  Output file name.
  * \param genome    Genome containing ref sequences (IDs and lengths are reported in VCF header).
  * \param vars      Variants to be output.
  * \param label     Header for the genotype column.
  * \returns         Number of exported variants.
  */
unsigned 
writeVcf (
  const std::string filename,
  const GenomeReference& genome,
  const std::vector<Variant>& vars,
  const std::string label
);

/** Read mutation map (clone x mutation) from a CSV file. */
int readMutMapFromCSV(
  std::map<std::string, std::vector<bool>> &mm,
  const std::string &filename
);

/** DEPRECATED! Functionality shifted to vario::VariantStore::generateGermlineVariants
    Generate variant loci in a given genome based on evolutionary model.
    Nucleotide substitution probabilities guide selection of loci. */
// std::vector<Variant> generateGermlineVariants(
//   const int num_variants,
//   const GenomeReference& genome,
//   GermlineSubstitutionModel& model,
//   RandomNumberGenerator<>&,
//   const bool infinite_sites = false
// );

/** Generate variant loci in a given genome based on evolutionary model.
    Loci are selected randomly. */
std::vector<Variant> generateVariantsRandomPos(
  const int num_variants,
  const GenomeReference& genome,
  GermlineSubstitutionModel& model,
  RandomNumberGenerator&,
  const bool infinite_sites = false
);
/** Apply variants to a given reference sequence */
void applyVariants(
  GenomeReference&,
  const std::vector<Variant>&);
/** Apply variants to a given reference sequence */
void applyVariants(
  GenomeReference&,
  const std::vector<Variant>&,
  const std::vector<Genotype>&);

// TODO: is this functionality still useful?
/** Apply variants to a given reference sequence. streaming modified genome to output. */
// void applyVariantsStream(
//   const Genome &ref_genome,
//   const std::vector<Mutation> &mutations,
//   const std::vector<Variant> &variants,
//   std::ostream &outstream,
//   short len_line = 60);

/** Print CopyNumberVariant details (BED format). */
std::ostream& operator<<(std::ostream& lhs, const CopyNumberVariant& cnv);


} /* namespace vario */

#endif /* VARIO_H */
