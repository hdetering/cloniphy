#ifndef BULKSAMPLEGENERATOR_H
#define BULKSAMPLEGENERATOR_H

#include "../bamio.hpp"

namespace bamio {

/* global typedefs to deal with SeqAn-internal types */
typedef seqan::StringSet<seqan::CharString> TNameStore;
typedef seqan::NameStoreCache<TNameStore> TNameStoreCache;
typedef seqan::BamIOContext<TNameStore, TNameStoreCache, seqan::Dependent<> > TBamContext;

/** Generates multi-sample NGS data */
class BulkSampleGenerator {

private:
  /** Flag that indicates if reference sequences have been set. */
  bool has_refseqs;
  /** Flag that indicates if clone genomes have been exported. */
  bool has_clone_genomes;
  /** The set of reference sequences to be included in output header. */
  std::map<std::string, seqio::TCoord> m_map_ref_len;
  /** Index of genomic segments for each clone and chromosome. */
  std::map<std::string, std::map<std::string, seqio::TSegMap>> m_map_clone_chr_seg;
  /** Total lengths of clone genomes (sum of SegmentCopies + padding) */
  std::map<std::string, unsigned long long> m_map_clone_len;
  /** Filenames of reference FASTA files, along with seq length, for each clone. */
  std::map<boost::filesystem::path, unsigned long long> m_map_fasta_len;
  /** Stores for each FASTA file the number of contained sequences. */
  std::map<boost::filesystem::path, unsigned> m_map_fasta_nseq;

public:
  /** Default c'tor. */
  BulkSampleGenerator();

  /** Extract the context of reference sequences (name, length) from reference genome. */
  void
  initializeRefSeqs (
    const seqio::GenomeReference& ref_genome
  );

  /** Attach reference sequences as BAM context to existing BAM header.
    * Requires reference sequences to be set (via initializeRefSeqs()).
    *
    * \param out_header  output param, BamHeader to be initialized.
    * \param out_context output param, BamIOContext to be initialized.
    * \returns true if successful, false on error.
    */
  bool
  generateBamHeader (
    seqan::BamHeader& out_header,
    TBamContext& out_context
  );

  /** Generate a set of bulk sequencing samples.
    *
    * \param mtx_sample_clone  Sampling matrix, contains clone abundances for each sample.
    * \param var_store         VariantStore containing germline and somatic variants.
    * \param vec_clone_lbl     Clone labels, corresponding to positions in sampling matrix.
    * \param path_fasta        Directory with clone genome sequences (tiled by copy number state).
    * \param path_bam          Directory to which reads (BAM) will be output.
    * \param seq_read_len      Read length (passed on to read simulator).
    * \param seq_frag_len_mean Insert size mean (passed on to read simulator).
    * \param seq_read_len_sd   Insert size std dev (passed on to read simulator).
    * \param seq_coverage      Haploid total sequencing coverage.
    * \param art_bin           Binary of read simulator to be called.
    * \param rng               Random number generator.
    */
  void
  generateBulkSamples (
    const std::map<std::string, std::map<std::string, double>> mtx_sample_clone,
    const vario::VariantStore var_store,
    const boost::filesystem::path path_fasta,
    const boost::filesystem::path path_bam,
    const unsigned seq_read_len,
    const unsigned seq_frag_len_mean,
    const unsigned seq_frag_len_sd,
    const double seq_coverage,
    const std::string art_bin,
    RandomNumberGenerator<>& rng
  );

  /** Generate sequencing reads for a bulk seq sample. 
    * Output files naming convention: 
    *   <bulk_sample>.<clone>.<copy_number>.sam
    *
    * \param path_fasta       Directory with clone genome sequences (tiled by copy number state).
    * \param path_bam         Directory to which reads (BAM) will be output.
    * \param lbl_sample       Label of the bulk sample (bulk_sample above).
    * \param map_clone_weight Clone weights, used to calculate read coverage for each clone.
    * \param seq_coverage     Haploid total sequencing coverage.
    * \param art              ArtWrapper object, called to generate reads from genomes.
    */
  void
  generateBulkSeqReads (
    const boost::filesystem::path path_fasta,
    const boost::filesystem::path path_bam,
    const std::string lbl_sample,
    const std::map<std::string, double> map_clone_weight,
    const double seq_coverage,
    ArtWrapper& art
  );

  /** Merge generated reads to form a bulk sample.
    * Output files naming convention: 
    *   <bulk_sample>.sam
    *
    * \param path_bam   directory containg read BAMs, output will be written there
    * \param lbl_sample label of the bulk sample (bulk_sample above)
    * \param vec_rg     read groups to include in the output BAM header
    * \param var_store  VariantStore containing germline and somatic variants.
    * \param rng        Random number generator.
    * \returns          true on success, false on error
    */
  bool
  mergeBulkSeqReads (
    const boost::filesystem::path path_bam,
    const std::string lbl_sample,
    const std::vector<seqan::BamHeaderRecord>& vec_rg,
    const vario::VariantStore var_store,
    RandomNumberGenerator<>& rng
  );

  /** Handle a tiled BAM file (representing a genomic region with consistent copy number). 
    * Reference sequence IDs are expected to follow the naming convention:
    *   <chromosome>_<start>_<end>_<padding>
    * Carries out the following tasks:
    * - transform mapping coordinates (local to global)
    * - discard reads in padded sequence
    * - spike-in variants
    *
    * \param bam_out   transformed reads are written to this BAM output file
    * \param bam_in    BAM file containing input reads (local to genomic tiles)
    * \param id_rg     read group ID to associate reads with
    * \param var_store VariantStore containing variants for spike-in.
    * \param rng       Random number generator.
    * \returns         true on success, false on error
    */
  bool
  transformBamTile (
    seqan::BamFileOut& bam_out,
    seqan::BamFileIn& bam_in,
    const std::string id_clone,
    const vario::VariantStore& var_store,
    RandomNumberGenerator<>& rng
  );

  /** Spike in mutations into reads.
    * 1. Select genomic SegmentCopy from interval map.
    * 2. Get variants associated with SegmentCopy from VariantStore.
    * 3. Identify variants overlapping read pair and apply them.
    *  
    * \param read1      First read in pair.
    * \param read2      Second read in pair.
    * \param segments   Interval map of genomic SegmentCopies.
    * \param var_store  VariantStore providing Variants-to-SegmentCopy mapping.
    * \param selector   Random selector (pick random collection element).
    * \returns          True on success, false on error.
    */
  bool
  mutateReadPair (
    seqan::BamAlignmentRecord& read1, 
    seqan::BamAlignmentRecord& read2,
    const seqio::TSegMap& segments,
    const vario::VariantStore& var_store,
    random_selector<>& selector
  );

  /** Add read groups to a given BamHeader. */
  void
  addReadGroups (
    seqan::BamHeader& header, 
    const std::vector<seqan::BamHeaderRecord>& vec_rg
  );

  /** Generate read groups to identify clones in bulk samples. */
  void
  generateReadGroups (
    std::vector<seqan::BamHeaderRecord>& vec_rg,
    const std::string id_sample, 
    const std::vector<std::string>& vec_tag_id,
    const std::string tag_pl,
    const std::string tag_pu
  );

  /** Exports genomes for individual clones, respecting individual CN states. 
    * \param map_lbl_gi     GenomeInstances to be exported, indexed by labels.
    * \param reference      Reference genome (contains actual DNA sequences).
    * \param padding        Genomic fragments will be padded by this number of 'A's.
    * \param min_len        Genomic fragments shorter than this will not be exported. 
    * \param path_fasta     Output folder for FASTA files.
    * \param path_bed       Output folder for BED files.
    * \param map_clone_glen Output param: total genome lengths (sum of all SegmentCopies) for each clone.
    */
  void
  writeCloneGenomes (
    const std::map<std::string, seqio::GenomeInstance> map_lbl_gi,
    const seqio::GenomeReference reference,
    unsigned padding,
    unsigned min_len,
    const boost::filesystem::path path_fasta,
    const boost::filesystem::path path_bed
  );

  /** Writes expected absolute copy number (CN) states for each bulk sample to a BED file.
    * Output: One BED file for each sample. 
    *         Naming convention: <sample>.cn.bed
    *
    * \param   mtx_sample sampling matrix, contains clone abundances for each sample.
    * \param   map_lbl_gi GenomeInstances making up bulk samples, indexed by clone labels.
    * \param   path_out   Output directory for BED files.
    * \returns true on success, false on error
    */
  bool
  writeBulkCopyNumber (
    const std::map<std::string, std::map<std::string, double>> mtx_sample,
    const std::map<std::string, seqio::GenomeInstance> map_lbl_gi,
    const boost::filesystem::path path_out
  ) const;

  /** Write genomic sequence of a GenomeInstance to FASTA file(s). 
    * Output:
    * - one FASTA file for each copy number (CN) state.
    * - single BED file with CN state by genomic region
    *
    * \param genome     GenomeInstance for which to export sequence.
    * \param reference  GenomeReference containing DNA sequence.
    * \param label      Filename prefix for output files. Will append "_<CN>.fa".
    * \param padding    Genomic fragments will be padded by this number of 'A's on both ends.
    * \param min_len    Genomic fragments shorter than this will not be exported. 
    * \param path_fasta Path to output FASTA file to.
    * \param path_bed   Path to output BED file to.
    */
  void
  writeFastaTiled (
    const seqio::GenomeInstance genome,
    const seqio::GenomeReference reference,
    const std::string label,
    unsigned int padding,
    unsigned int min_len,
    const boost::filesystem::path path_fasta,
    const boost::filesystem::path path_bed
  );
};

}

#endif /* BULKSAMPLEGENERATOR_H */