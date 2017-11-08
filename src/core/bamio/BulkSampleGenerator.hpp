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
  /** Allele counts of SNVs indexed by clone, SNV id. */
  std::map<std::string, std::map<int, vario::VariantAlleleCount>> m_map_clone_snv_vac;
  /** Allele frequencies of SNVs indexed by SNV id. */
  std::map<int, double> m_map_snv_vaf;

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
    * \param path_bam          Directory to which allele counts (BED) will be output.
    * \param seq_coverage      Total sequencing coverage (haploid).
    * \param seq_error         Sequencing error (per base). Only applied when generating read counts.
    * \param seq_read_gen      Generate reads? (false: generate read counts)
    * \param seq_use_vaf       Spike in variants according to VAFs? (breaks haplotypes!)
    * \param seq_read_len      Read length (passed on to read simulator).
    * \param seq_frag_len_mean Insert size mean (passed on to read simulator).
    * \param seq_read_len_sd   Insert size std dev (passed on to read simulator).
    * \param art_bin           Binary of read simulator to be called.
    * \param rng               Random number generator.
    */
  void
  generateBulkSamples (
    const std::map<std::string, std::map<std::string, double>> mtx_sample_clone,
    const vario::VariantStore var_store,
    const boost::filesystem::path path_fasta,
    const boost::filesystem::path path_bam,
    const boost::filesystem::path path_bed,
    const double seq_coverage,
    const double seq_error,
    const bool seq_read_gen,
    const bool seq_use_vaf,
    const unsigned seq_read_len,
    const unsigned seq_frag_len_mean,
    const unsigned seq_frag_len_sd,
    const std::string art_bin,
    RandomNumberGenerator<>& rng
  );

  /**
   * Generate read counts for variant loci for each sample.
   * 
   * Output files naming convention:
   *   <sample>.vars.csv
   * 
   * Output file column format:
   *   <mutId>,<total_reads>,<alternative_reads>
   * 
   * \param path_out      Path at which to write output files.
   * \param lbl_sample    Label for the sample.
   * \param seq_coverage  Sequencing depth.
   * \param seq_error     Sequencing error (per base).
   * \returns             true on success, false on error
   */
  bool
  generateReadCounts (
    const boost::filesystem::path path_out,
    const std::string lbl_sample,
    const double seq_coverage,
    const double seq_error,
    RandomNumberGenerator<>& rng
  ); 

  /** 
   * Generate sequencing reads for a bulk seq sample. 
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
    * \param path_bam     Directory containg read BAMs, output will be written there.
    * \param lbl_sample   Label of the bulk sample (bulk_sample above).
    * \param vec_rg       Read groups to include in the output BAM header.
    * \param var_store    VariantStore containing germline and somatic variants.
    * \param seq_use_vaf  Spike in variants according to VAFs? (breaks haplotypes!)
    * \param rng          Random number generator.
    * \returns            true on success, false on error
    */
  bool
  mergeBulkSeqReads (
    const boost::filesystem::path path_bam,
    const std::string lbl_sample,
    const std::vector<seqan::BamHeaderRecord>& vec_rg,
    const vario::VariantStore var_store,
    const bool seq_use_vaf,
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
    * \param bam_out      transformed reads are written to this BAM output file
    * \param bam_in       BAM file containing input reads (local to genomic tiles)
    * \param id_rg        read group ID to associate reads with
    * \param var_store    VariantStore containing variants for spike-in.
    * \param rng          Random number generator.
    * \param map_var_cvg  Output: total read count for each variant.
    * \param map_var_alt  Output: alternative allele read count for each variant.
    * \returns            true on success, false on error
    */
  bool
  transformBamTileSeg (
    seqan::BamFileOut& bam_out,
    seqan::BamFileIn& bam_in,
    const std::string id_clone,
    const vario::VariantStore& var_store,
    RandomNumberGenerator<>& rng,
    std::map<std::string, unsigned>& map_var_cvg,
    std::map<std::string, unsigned>& map_var_alt
  );

  /** Using mutateReadPairOpt code, but directly (without function call). */
  bool
  transformBamTileVaf (
    seqan::BamFileOut& bam_out,
    seqan::BamFileIn& bam_in,
    const std::string id_clone,
    const vario::VariantStore& var_store,
    RandomNumberGenerator<>& rng,
    std::map<std::string, unsigned>& map_var_cvg,
    std::map<std::string, unsigned>& map_var_alt
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
  mutateReadPairSeg (
    seqan::BamAlignmentRecord& read1, 
    seqan::BamAlignmentRecord& read2,
    const seqio::TSegMap& segments,
    const vario::VariantStore& var_store,
    random_selector<>& selector
  );

  /** Spike in mutations into reads.
    * 1. Identify first variant downstream of read pair start.
    * 2. Apply variants overlapping read pair according to their expected VAF.
    *  
    * \param read1       First read in pair.
    * \param read2       Second read in pair.
    * \param map_pos_snv SNV ids indexed by position (for quick lookup).
    * \param map_id_snv  SNV Variant objects indexed by id.
    * \param r_dbl       Random function (value decides if read receives mutation).
    * \param r1_begin    Start coordinate of read1.
    * \param r2_begin    Start coordinate of read2.
    * \param r1_end      End coordinate of read1.
    * \param r2_end      End coordinate of read2.
    * \returns           True on success, false on error.
    */
  bool
  mutateReadPairVaf (
    seqan::BamAlignmentRecord& read1, 
    seqan::BamAlignmentRecord& read2,
    const std::map<seqio::TCoord, std::vector<int>>& map_pos_snv,
    const std::map<int, vario::Variant>& map_id_snv,
    std::function<double()>& r_dbl,
    const int r1_begin,
    const int r2_begin,
    const int r1_end,
    const int r2_end
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

  /** 
   * Calculate allele counts at variant positions for all clones.
   * Initializes m_map_clone_snv_vac, m_map_snv_vaf.
   *
   * \param map_clone_ccf  Cancer cell fraction (CCF) for each clone in the sample.
   * \param var_store      Variant store keeping track of SNV -> segment copy mappings.
   * \returns              True on success, false on error.
   */
  bool
  initAlleleCounts (
    const std::map<std::string, double> map_clone_ccf,
    const vario::VariantStore& var_store
    //vario::TMapChrPosVaf& out_vars
  );

  /** 
   * Write expected variant allele frequency (VAF) for a sample to CSV file.
   * Expected read counts for each variant are calculated as:
   *   (\sum_i A_i / R_i * F_i) * C
   * for each clone i, where
   *   A_i := number of alternative allele copies
   *   R_i := number of reference allele copies
   *   F_i := frequency of clone i in the sample
   *   C   := depth of coverage for the sample
   * 
   * Output format:
   *   id_snv,chr,pos,alt
   * 
   * \param ofs_bed_out  Output stream to write output to.
   * \param cvg_depth    Coverage depth during read simulation.
   * \param var_store    Variant store, contains details about variants.
   * \returns            Number of output variants on success, -1 on error.
   */
  int
  writeExpectedReadCounts (
    std::ofstream& ofs_bed_out,
    const int cvg_depth,
    const vario::VariantStore& var_store
  ) const;
};

}

#endif /* BULKSAMPLEGENERATOR_H */