#ifndef SEQIO_H
#define SEQIO_H

#include "random.hpp"
#include "stringio.hpp"
#include <algorithm>
#include <boost/circular_buffer.hpp>
#include <boost/format.hpp>
#include <boost/icl/interval.hpp>
#include <boost/icl/interval_map.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <cassert>
#include <fstream>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <memory> // unique_ptr, shared_ptr, weak_ptr
#include <sstream>
#include <string>
#include <vector>

/** Handles sequence files (read, write, index). */
namespace seqio {

/** DEPRECATED! Use TCoord instead. */
typedef unsigned long ulong;
/** Represents genomic coordinates. */
typedef unsigned long TCoord;
/** Represents genomic regions (chromosome, start, end). */
typedef std::tuple<std::string, ulong, ulong> TRegion;
/** Set of valid nucleotides. */
enum Nuc { A, C, G, T, N };

/**
 * Returns a canonical index (corresponding to ::Nuc) for a nucleotide.
 *
 * \returns index position of nucleotide, -1 for unknown characters
 */
inline short nuc2idx (char nuc) {
  switch (nuc) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
  }
  return -1; // default for unknown chars
}

/**
 * Map nucleotide char to a canonical index (corresponding to ::Nuc).
 *
 * \returns Nucleotide char for defined indices, `?` otherwise
 */
inline char idx2nuc (short idx) {
  switch (idx) {
    case 0: return 'A';
    case 1: return 'C';
    case 2: return 'G';
    case 3: return 'T';
    case 4: return 'N';
  }
  return '?'; // default for unknown chars
}

/**
 * Get reverse complement of a nucleotide.
 */
inline char rev_comp (char nuc) {
   switch (nuc) {
     case 'A': case 'a': return 'T';
     case 'C': case 'c': return 'G';
     case 'G': case 'g': return 'C';
     case 'T': case 't': return 'A';
   }
   assert(false);
   return '?'; // default for unknown chars
}

/**
 * Get reverse complement of a DNA sequence.
 */
inline std::string rev_comp (std::string dna) {
  std::transform(dna.begin(), dna.end(), dna.begin(),
                 [](char c) { return rev_comp(c); });
  return dna;
}

struct SeqRecord
{
  std::string id;          /** identifier */
  std::string description; /** sequence description (everything after first space in ID line) */
  std::string seq;         /** actual sequence */
  std::string id_ref;      /** identifier in reference genome (ploidy) */
  short chr_copy;          /** chromosome copy (0 for haploid) */
  SeqRecord(const std::string, const std::string, const std::string&);
  ~SeqRecord();
};

/** Represents a genomic location */
struct Locus
{
  ulong       idx_record; // index of genomic 
  TCoord      start;      // global start position in sequence (0-based, inclusive)
  TCoord      end;        // global end position in sequence (0-based, exclusive)
  std::string id_ref;     // seq id in ref genome

  Locus();
  Locus(std::string id, TCoord start, TCoord end);
  ~Locus();
};

/** Represents a copy of a reference sequence segment.
 *
 *  Inititally, each allele (n=2 in the normal, diploid case) is a copy of
 *  the reference sequence. Further segment copies arise from CNV events.
 */
struct SegmentCopy {
  /** Universally-unique identifier to globally refer to the object. */
  boost::uuids::uuid id;
  /** Start position (0-based, inclusive) in reference chromosome */
  TCoord ref_start;
  /** End position (0-based, exclusive) in reference chromosome */
  TCoord ref_end;
  /** Germline source allele (maternal: A, paternal: B). */
  char gl_allele;

  /** default c'tor */
  SegmentCopy();
  /** dummy c'tor to avoid UUID init (to avoid performance overhead) */
  SegmentCopy(int);
  /** default d'tor */
  ~SegmentCopy();

  /** Create SegmentCopy for given coordinates.
   *  \param start      Reference bp position where segment copy begins.
   *  \param end        Reference bp position where segment copy ends.
   *  \param gl_allele  Germline source allele ('A','B') from which segment copy originates.
   */
  SegmentCopy (
    const TCoord start, 
    const TCoord end,
    const char gl_allele
  );
};

/** Compare two segment copies (equality).
 *  If IDs match, segment copies are identical.
 */
bool operator==(const SegmentCopy& lhs, const SegmentCopy& rhs);

/** Compare two segment copies (order).
 *  Order by UUID (necessary to avoid collisions in boost::icl::set).
 */
bool operator<(const SegmentCopy& lhs, const SegmentCopy& rhs);

/** Stores SegmentCopies overlapping genomic intervals. */
typedef boost::icl::interval_map<TCoord, std::set<SegmentCopy>> TSegMap;
/** Stores SegmentCopies as interval set */
typedef boost::icl::interval_set<TCoord> TSegSet;

/** Represents the reference genome version of a chromosome.
 *
 *  Structurally a ChromosomeReference need not be a contiguous sequence,
 *  but may be composed of shorter elements (e.g., exons, targets).
 */
struct ChromosomeReference {
   /** chromosome identifier in reference genome. */
   std::string id;
   /** total length */
   unsigned long length;
   /** sequence records belonging to this chromosome, indexed by start position */
   std::map<ulong, std::shared_ptr<SeqRecord>> map_start_rec;

   /** default c'tor */
   ChromosomeReference();

   /** Get nucleotide at position */
   char getNucAt(unsigned long pos);
 };

/** Documents a modification introduced to a SegmentCopy.
 *  Structure: (new, old, start, end)
 *   - new: UUID of newly created SegmentCopy
 *   - old: UUID of existing SegmentCopy to be replaced
 *   - start: start coordinate within old SegmentCopy covered by the new one
 *   - end: end coordinate within old SegmentCopy covered by the new one
 */
typedef std::tuple<boost::uuids::uuid, boost::uuids::uuid, unsigned long, unsigned long> seg_mod_t;

/** Represents a physical representation of a chromosome.
 *
 *  Structurally a ChromosomeInstance is a chain of SegmentCopies,
 *  which represent regions of the reference sequence.
 */
struct ChromosomeInstance {
  /** Real length (in bp) of ChromosomeInstance */
  unsigned long length;
  /** SegmentCopies that are associated with this ChromosomeInstance */
  std::list<SegmentCopy> lst_segments;
  
  /** default c'tor */
  ChromosomeInstance();

  /** Generate ChromsomeInstance from reference chromosome.
   *  A single SegmentCopy will be created comprising the whole chromosome.
   *  \param chr_ref    Reference chromosome of which to create an instance.
   *  \param gl_allele  Germline source allele, will be assigned to segment copies.
   */
  ChromosomeInstance (
    const ChromosomeReference chr_ref,
    const char gl_allele
  );

  /** Identify SegmentCopies overlapping a given locus. */
  std::vector<SegmentCopy> 
  getSegmentCopiesAt (
    TCoord ref_pos
  );

  /** Amplify a region of this ChromosomeInstance by creating new SegmentCopies
   *  \param start_rel    Relative start coordinate of deletion (fraction of chromosome length).
   *  \param len_rel      Relative length of region to delete (fraction of chromosome length).
   *  \param is_forward   true: Deletion affects region towards 3' end from start_rel; false: towards 5' end.
   *  \param is_telomeric true: Deletion includes chromosome end (facilitates end coordinate calculation)
   *  \param is_deletion  true: Deletion event; false: Amplification event
   */
  std::vector<seg_mod_t> 
  amplifyRegion (
    double start_rel,
    double len_rel,
    bool is_forward,
    bool is_telomeric
  );

  /** Delete a region of this ChromosomeInstance by creating new SegmentCopies.
   *  \param start_rel    Relative start coordinate of deletion (fraction of chromosome length).
   *  \param len_rel      Relative length of region to delete (fraction of chromosome length).
   *  \param is_forward   true: Deletion affects region towards 3' end from start_rel; false: towards 5' end.
   *  \param is_telomeric true: Deletion includes chromosome end (facilitates end coordinate calculation)
   */
  std::vector<seg_mod_t> 
  deleteRegion (
    double start_rel,
    double len_rel,
    bool is_forward,
    bool is_telomeric
  );

  /** Create a copy of an existing ChromosomeInstance.
   *
   *  Copies each SegmentCopy comprising the existing ChromosomeInstance.
   *  \param ci_old Existing ChromosomeInstance to be copied
   *  \param out_seg_mods Output parameter: Tuples of modifications to SegmentCopies (id_new, id_old, start_old, end_old)
   */
  void 
  copy (
    std::shared_ptr<ChromosomeInstance> ci_old, 
    std::vector<seg_mod_t>& out_seg_mods
  );

  /** Build interval map of SegmentCopies.
   *  Used to answer questions like: 
   *    "Which SegmentCopies overlap interval 'chr1:10,000-10,500'?"
   *  \param imap_segments output param: interval map of SegmentCopies
   *  \returns true on success, false on error
   */
   bool
   indexSegmentCopies (
     TSegMap& imap_segments
  ) const;
};

/** Stores a set of SeqRecords along with indexing information. */
struct GenomeReference
{
  unsigned num_records;
  unsigned length;                        /** total length of all sequences */
  unsigned masked_length;                 /** length of unmasked (non-'N') positions */
  // TODO: remove ploidy (modeled on SegmentCopy level)
  short ploidy;                           /** number of copies for each chromosome */
  // TODO: move SeqRecords below ChromosomeReference level
  std::vector<std::shared_ptr<SeqRecord>> records;
  /** Reference chromosomes that make up the genome */
  std::map<std::string, std::shared_ptr<ChromosomeReference>> chromosomes;
  /** Vector of chromosome IDs (paired with length vector) */
  std::vector<std::string> vec_chr_id;
  /** Vector of chromosome lengths (paired with chromosome ID vector) */
  std::vector<TCoord> vec_chr_len;

  std::vector<unsigned> vec_start_chr;    /** cumulative start positions of sequences */
  std::vector<unsigned> vec_start_masked; /** cumulative start positions of unmasked regions */
  std::vector<unsigned> vec_cumlen_masked; /** cumulative lengths of unmasked regions */
  double nuc_freq[4];                     /** nucleotide frequencies */
  /** absolute bp positions indexed by nucleotide */
  //std::map<char, std::vector<long> > map_nuc_pos;
  std::vector<std::vector<long> > nuc_pos;
  /** absolute bp positions indexed by tri-nucleotides */
  std::map<std::string, std::vector<long> > map_3mer_pos;

  /** default c'tor */
  GenomeReference();
  /** load reference sequences from FASTA file */
  GenomeReference(const char* fn_fasta);

  /** Adds a ChromosomeReference to this GenomeReference. */
  void addChromosome(std::shared_ptr<ChromosomeReference> sp_chr);

  /** Simulate DNA seq of given length and nuc freqs. */
  void generate (
    const unsigned long,
    const std::vector<double>,
    RandomNumberGenerator&
  );

  /**
   * Generate random reference genome based on nucleotide frequencies.
   *
   * \param total genome length including all sequences
   * \param number of chromosomes/sequences to generate
   * \param nucleotide frequencies (A,C,G,T)
   * \param object to generate random numbers
   */
  void generate (
    const unsigned long,
    const unsigned short,
    const std::vector<double>,
    RandomNumberGenerator&
  );

  /**
    * Generate random genome with given number of fragments, mean len, sd len, nuc freqs.
    *
    * \param num_frags  number of chromosomes/sequences to generate
    * \param mean_len   mean sequence length
    * \param sd_len     standard deviation of sequence length
    * \param nuc_freqs  nucleotide frequencies (A,C,G,T)
    * \param rng object to generate random numbers
    */
  void generate (
    const unsigned num_frags,
    const unsigned long mean_len,
    const unsigned long sd_len,
    const std::vector<double> nuc_freqs,
    RandomNumberGenerator& rng
  );

  /**
   * Index genome for easier access to sequences and individual nucleotides.
   *
   * Indices are generated to translate between global and local coordinates
   * as well as to quickly locate
   *  - unmasked (non-"N") regions
   *  - positions having a given nucleotide
   *  - positions having a given tri-nucleotide
   */
  void indexRecords();

  /**
   * Get chromosome and local position for global position
   */
  Locus getLocusByGlobalPos(long) const;

  /**
   * Get absolute coordinates for a relative position in unmasked part of the genome
   */
  Locus getAbsoluteLocusMasked(double) const;

  /** 
   * Get DNA sequence for a genomic region.
   * 
   * In the case of discontinuous reference (e.g. exome, targeted seq),
   * the output may be shorter than coordinates indicate.
   * 
   * \param id_chr  reference chromosome id
   * \param start   start coordinate
   * \param end     end coordinate
   * \param seqs    output parameter, sequences located within given coordinates
   */
   void 
   getSequence (
     const std::string id_chr,
     const TCoord start,
     const TCoord end,
     std::map<TCoord, std::string>& seqs
   ) const;

   /** 
    * Get nucleotide sequence at given locus.
    *
    * \param id_chr     reference chromosome id
    * \param pos_start  start coordinate (inclusive)
    * \param pos_end    end coordinate (exclusive)
    * \param seq        output parameter, sequence within given coordinates
    * \returns          true on success, false on error
    */
   bool
   getSequence (
     const std::string id_chr,
     const TCoord pos_start,
     const TCoord pos_end,
     std::string& seq
   ) const;
};

/** Reads sequences from file. */
void readFasta(const char*, std::vector<std::shared_ptr<SeqRecord>>&);
/** Reads sequences from istream. */
void readFasta(std::istream&, std::vector<std::shared_ptr<SeqRecord>>&);
/** Writes sequences to file. */
int writeFasta(const std::vector<std::shared_ptr<SeqRecord>>&, const std::string fn, int len_line = 60);
/** Writes sequences to ostream. */
int writeFasta(const std::vector<std::shared_ptr<SeqRecord>>&, std::ostream& os, int len_line = 60);
/** Generate an index for a FASTA file containing multiple sequences */
void indexFasta(const char*);
/** Generate random DNA sequence */
unsigned long generateRandomDnaSeq(
  std::string &seq,
  const unsigned long total_len,
  const std::vector<double> nuc_freqs,
  RandomNumberGenerator &rng);
/** Simulate allelic dropout events, masking parts of genome as 'N's. */
void simulateADO_old(const std::string, const float, const int, std::function<double()>&);
/** Simulate allelic dropout events, masking parts of genome as 'N's. */
void simulateADO(const std::string, const unsigned, const float, const int, std::function<double()>&);
/** Convert a nucleotide char into Nuc */
Nuc charToNuc(const char);
/** Convert Nuc into a nucleotide char */
char nucToChar(const Nuc);
/** Convert Nuc into a nucleotide string */
std::string nucToString(const Nuc);
/** Mutate a nucleotide by an offset. */
char shiftNucleotide(const char, const int);

/** Utility class to encapsulate allele-specific copy number state as two double values. */
struct AlleleSpecCopyNum {
  /** Copy number for allele A. */
  double count_A;
  /** Copy number for allele B. */
  double count_B;

  /** default c'tor */
  AlleleSpecCopyNum() : count_A(0.0), count_B(0.0) {}

  /** Add two copy number states. */
  AlleleSpecCopyNum& operator+=(const AlleleSpecCopyNum& rhs) {
    this->count_A += rhs.count_A;
    this->count_B += rhs.count_B;
    return *this;
  }
};
/** Compare two copy number states. */
bool operator==(const AlleleSpecCopyNum& lhs, const AlleleSpecCopyNum& rhs);

/** Represents the incarnation of a reference genome.
 *
 *  A GenomeInstance consists of a set of ChromosomeInstances,
 *  which in turn consist of SegmentCopies.
 */
struct GenomeInstance {
  /** Stores the ChromosomeInstances that make up this GenomeInstance. */
  std::vector<std::shared_ptr<ChromosomeInstance>> vec_chr;
  /** Stores for each ChromosomeInstance its real length. Used to randomly select a chromosome by length. */
  std::vector<unsigned long> vec_chr_len;
  /** ChromosomeInstances indexed by reference id. */
  std::map<std::string, std::vector<std::shared_ptr<ChromosomeInstance>>> map_id_chr;
  /** SegmentCopy interval maps indexed by chromosome id. */
  std::map<std::string, TSegMap> map_chr_seg;

  /** default c'tor */
  GenomeInstance();

  /** Create GenomeInstance from GenomeReference object
   *
   *  Create a diploid genome by initializing 2 ChromosomeInstances
   *  for each ChromosomeReference.
   */
  GenomeInstance(const GenomeReference& g_ref);

  /**
   * Create GenomeInstance as copy of existing object.
   * \param g_inst          Existing GenomeInstance to be copied.
   * \param out_vec_seg_mod Output param: mapping from original to new segment copy IDs.
   */
  GenomeInstance(
    const GenomeInstance& g_inst,
    std::vector<seg_mod_t>& out_vec_seg_mod
  );

  /** Adds a ChromosomeInstance to this GenomeInstance.
   *  \param sp_chr shared pointer to ChromosomeInstance
   *  \param id_chr reference ID under which to index ChromsomeInstance
   */
  void
  addChromosome(
    std::shared_ptr<ChromosomeInstance> sp_chr,
    std::string id_chr
  );

  /** Remove a ChromosomeInstance from this GenomeInstance.
   *  \param sp_chr shared pointer to ChromosomeInstance
   *  \param id_chr reference ID under which to index ChromsomeInstance
   */
  void
  deleteChromosome(
    std::shared_ptr<ChromosomeInstance> sp_chr,
    std::string id_chr
  );

  /**
   * Perform Whole Genome Duplication.
   *
   * Duplication is carried out by duplicating all ChromosomeInstances.
   * \param out_vec_seg_mod Output parameter: vector of SegmentCopy modifications (id_new, id_old, start_old, end_old)
   */
  void
  duplicate(std::vector<seg_mod_t>& out_vec_seg_mod);

  /** Identify SequenceCopies overlapping a given locus. */
  std::vector<SegmentCopy>
  getSegmentCopiesAt (
    std::string id_chromosome,
    unsigned long ref_pos
  );

  /**DEPRECATED!
   * Map genomic regions to their absolute copy number states.
   *  \param out_map_cn_loci map storing for each copy number state a collection of genomic loci.
   */
  void
  getCopyNumberStates (
    std::map<unsigned, std::vector<std::shared_ptr<Locus>>>& out_map_cn_loci
  );

  /** Map genomic regions for each chromosome to their absolute copy number states.
   *  \param map_chr_seg map storing for genomic segment its copy number state.
   *  \param scale       factor by which each segment increases CN (default: 1).
   */
  void
  getCopyNumberStateByChr (
    std::map<std::string, boost::icl::interval_map<TCoord, double>>& map_chr_seg,
    const double scale = 1.0
  ) const;

  /** Map genomic regions for each chromosome to their allele-specific copy number states.
   *  \param map_chr_seg map storing for genomic segment its copy number state.
   *  \param scale       factor by which each segment increases CN (default: 1).
   */
  void
  getCopyNumberStateByChr (
    std::map<std::string, boost::icl::interval_map<TCoord, AlleleSpecCopyNum>>& map_chr_seg,
    const double scale = 1.0
  ) const;

  /** Build interval map of SegmentCopies for each chromosome.
   *  Used to answer questions like: 
   *    "Which SegmentCopies overlap interval 'chr1:10,000-10,500'?"
   *  \returns true on success, false on error
   */
   bool
   indexSegmentCopies ();
};

/** Print GenomeInstance substructure (for debugging purposes). */
std::ostream& operator<<(std::ostream& lhs, const GenomeInstance& gi);
/** Print ChromosomeInstance substructure (for debugging purposes). */
std::ostream& operator<<(std::ostream& lhs, const ChromosomeInstance& ci);
/** Print SegmentCopy info (for debugging purposes). */
std::ostream& operator<<(std::ostream& lhs, const SegmentCopy& seg);

} /* namespace seqio */

#endif /* SEQIO_H */
