#ifndef GENOMEREFERENCE_H
#define GENOMEREFERENCE_H

#include "../random.hpp"
#include "ChromosomeReference.hpp"
#include "KmerProfile.hpp"
#include "Locus.hpp"
#include "SeqRecord.hpp"
#include <map>
#include <memory> // unique_ptr, shared_ptr, weak_ptr
#include <string>
#include <vector>

namespace seqio {

/** Stores a set of SeqRecords along with indexing information. */
struct GenomeReference
{
  unsigned num_records;
  unsigned length;                        /** total length of all sequences */
  unsigned masked_length;                 /** length of unmasked (non-'N') positions */
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
  /** load reference sequences from FASTA file, filter by Locus list. */
  GenomeReference(
    const char* fn_fasta,
    const std::map<std::string, std::vector<Locus>>& loci
  );

  /** Adds a ChromosomeReference to this GenomeReference. */
  void addChromosome(std::shared_ptr<ChromosomeReference> sp_chr);

  /** Simulate DNA seq of given length and nuc freqs. */
  void generate_nucfreqs (
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
  void generate_nucfreqs (
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
  void generate_nucfreqs (
    const unsigned num_frags,
    const unsigned long mean_len,
    const unsigned long sd_len,
    const std::vector<double> nuc_freqs,
    RandomNumberGenerator& rng
  );

  /**
    * Generate random genome with given number of fragments, mean len, sd len, nuc freqs.
    * The generated sequence will contain trinucleotides according to the specified
    * frequencies.
    *
    * \param num_frags   number of chromosomes/sequences to generate
    * \param mean_len    mean sequence length
    * \param sd_len      standard deviation of sequence length
    * \param kmer_prof   kmer frequencies
    * \param rng         object to generate random numbers
    */
  bool generate_kmer (
    const unsigned num_frags,
    const unsigned long mean_len,
    const unsigned long sd_len,
    const KmerProfile kmer_prof,
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

}; // struct GenomeReference

} // namespace seqio

#endif // GENOMEREFERENCE_H