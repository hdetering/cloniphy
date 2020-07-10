#include "../seqio.hpp"
#include "GenomeReference.hpp"
#include <cmath> // floor()

using namespace std;

namespace seqio {

GenomeReference::GenomeReference() : length(0), masked_length(0) {}

GenomeReference::GenomeReference(const char* filename)
: length(0), 
  masked_length(0)
{
  // read records from file
  readFasta(this->records, filename);
  // scan records for segmented sequence (naming convention: "CHR:START-END")
  for (auto & rec : this->records) {
    // TODO: does it make sense to expect genomic fragments in FASTA?
    // if (regex_search(rec.id, matches, rgx)) { // segment found
    //   string id_chr = matches[1];
    //   unsigned long start = atoi(matches[2]);
    //   unsigned long end = atoi(matches[3]);
    //
    //   // check if chromosome already exists
    //   auto kv = chromosomes.find(id_chr);
    //   if (kv == chromosomes.end()) { // entry not found
    //     ChromosomeReference chr_ref;
    //     chr_ref.id = rec.id;
    //     chromosomes[id_chr] =
    //   }
    //   else { // entry exists
    //   }
    // }
    // else { // treat sequence as chromosome
    // }
    shared_ptr<ChromosomeReference> p_chr_ref(new ChromosomeReference());
    p_chr_ref->id = rec->id;
    p_chr_ref->length = rec->seq.length();
    p_chr_ref->map_start_rec[0] = rec;
    // sanity check: chromosome IDs should be unique
    assert(this->chromosomes.count(p_chr_ref->id) == 0);
    this->addChromosome(p_chr_ref);
  }
  num_records = records.size();
}

GenomeReference::GenomeReference (
  const char* filename,
  const map<string, vector<Locus>>& map_chr_loci
)
: length(0), 
  masked_length(0) 
{
  // read records from file
  bool use_whitelist = ( map_chr_loci.size() > 0 );
  readFasta(this->records, filename, use_whitelist, map_chr_loci);
  // scan records
  for (auto & rec : this->records) {
    // if record ID is not present in loci list: skip record
    if (map_chr_loci.find(rec->id) == map_chr_loci.end())
      continue;
    // initialize new reference chromosome
    shared_ptr<ChromosomeReference> p_chr_ref(new ChromosomeReference());
    p_chr_ref->id = rec->id;
    p_chr_ref->length = rec->seq.length();
    for (const Locus loc : map_chr_loci.at(rec->id)) {
      string seq_id = stringio::format("%s:%lu-%lu", loc.id_ref, loc.start, loc.end);
      string seq = rec->seq.substr(loc.start, loc.end-loc.start);
      shared_ptr<SeqRecord> sp_rec(new SeqRecord(seq_id, "", seq));
      p_chr_ref->map_start_rec[loc.start] = sp_rec;
    }
    // sanity check: chromosome IDs should be unique
    assert(this->chromosomes.count(p_chr_ref->id) == 0);
    this->addChromosome(p_chr_ref);
  }
  num_records = records.size();
}

GenomeReference::~GenomeReference () {}

void GenomeReference::addChromosome(shared_ptr<ChromosomeReference> sp_chr) {
  this->chromosomes[sp_chr->id] = sp_chr;
  this->vec_chr_id.push_back(sp_chr->id);
  this->vec_chr_len.push_back(sp_chr->length);
}

void GenomeReference::generate_nucfreqs (
  const unsigned long total_len,
  const vector<double> nuc_freqs,
  RandomNumberGenerator &rng
) {
  this->generate_nucfreqs(total_len, 1, nuc_freqs, rng);
}

/**
 * Generate random reference genome based on nucleotide frequencies.
 *
 * \param total_len total genome length including all sequences
 * \param num_chr number of chromosomes/sequences to generate
 * \param nuc_freqs nucleotide frequencies (A,C,G,T)
 * \param object to generate random numbers
 */
void GenomeReference::generate_nucfreqs (
  const unsigned long total_len,
  const unsigned short num_chr,
  const vector<double> nuc_freqs,
  RandomNumberGenerator &rng
) {
  // generate random genomic sequence
  string seq;
  unsigned long gen_len = 0;
  gen_len = generateRandomDnaSeq(seq, total_len, nuc_freqs, rng);

  // generate chromosome limits
  vector<unsigned long> chr_ends = { 0 };
  if (num_chr > 1) {
    unsigned long zero = 0;
    function<unsigned long()> rpos = rng.getRandomFunctionInt(zero, total_len);
    for (auto i=0; i<num_chr-1; ++i) {
      unsigned long p = rpos();
      chr_ends.push_back(p);
    }
    sort(chr_ends.begin(), chr_ends.end());
  }
  chr_ends.push_back(total_len);

  // split genome into chromosome sequences
  unsigned long cum_start = 0;
  string::iterator it_start = seq.begin();
  string::iterator it_end = seq.begin();
  for (auto i=0; i<num_chr; ++i) {
    it_start = it_end;
    advance(it_end, chr_ends[i+1]-chr_ends[i]);
    string id_chr = stringio::format("chr%d", i);
    shared_ptr<SeqRecord> sp_rec(new SeqRecord(id_chr, "random sequence", string(it_start, it_end)));
    this->records.push_back(sp_rec);
    // instantiate new referennce chromosome
    shared_ptr<ChromosomeReference> sp_chr(new ChromosomeReference());
    sp_chr->id = id_chr;
    sp_chr->length = it_end - it_start;
    sp_chr->map_start_rec[0] = sp_rec;
    this->addChromosome(sp_chr);
  }
  this->num_records = num_chr;
  this->length = gen_len;
}

/** generate random genome by given number of fragments, mean len, sd len, nuc freqs */
void GenomeReference::generate_nucfreqs (
  const unsigned num_seqs,
  const unsigned long mean_len,
  const unsigned long sd_len,
  const vector<double> nuc_freqs,
  RandomNumberGenerator& rng)
{
  // set sequence lengths
  vector<unsigned long> vec_seq_len(num_seqs);
  if (sd_len > 0) { // sample from Gamma distribution
    function<double()> rlen = rng.getRandomFunctionGammaMeanSd(mean_len, sd_len);
    generate(vec_seq_len.begin(), vec_seq_len.end(), rlen);
    sort(vec_seq_len.begin(), vec_seq_len.end(), std::greater<unsigned long>());
  }
  else { // fixed length for all sequences
    fill(vec_seq_len.begin(), vec_seq_len.end(), mean_len);
  }

  // simulate sequences
  unsigned idx_chr = 1;
  for (auto l : vec_seq_len) {
    string seq;
    generateRandomDnaSeq(seq, l, nuc_freqs, rng);
    string id_chr = stringio::format("chr%d", idx_chr++);
    shared_ptr<SeqRecord> sp_rec(new SeqRecord(id_chr, "random sequence", seq));
    sp_rec->id_ref = sp_rec->id;
    this->records.push_back(sp_rec);
    // instantiate new referennce chromosome
    shared_ptr<ChromosomeReference> sp_chr(new ChromosomeReference());
    sp_chr->id = id_chr;
    sp_chr->length = seq.length();
    sp_chr->map_start_rec[0] = sp_rec;
    this->addChromosome(sp_chr);
  }
  this->num_records = num_seqs;
}

/** 
 * Generate random genome by given
 * - number of fragments
 * - mean len
 * - sd len
 * - kmer freqs
 */
bool GenomeReference::generate_kmer (
  const unsigned num_seqs,
  const unsigned long mean_len,
  const unsigned long sd_len,
  const KmerProfile kmer_prof,
  RandomNumberGenerator& rng)
{
  // make sure KmerProfile has prefixes indexed
  assert( kmer_prof.has_idx_pfx );

  // GENERATE SEQUENCES
  //---------------------------------------------------------------------------

  // set sequence lengths
  vector<unsigned long> vec_seq_len(num_seqs);
  if (sd_len > 0) { // sample from Gamma distribution
    function<double()> rlen = rng.getRandomFunctionGammaMeanSd(mean_len, sd_len);
    generate(vec_seq_len.begin(), vec_seq_len.end(), rlen);
    sort(vec_seq_len.begin(), vec_seq_len.end(), std::greater<unsigned long>());
  }
  else { // fixed length for all sequences
    fill(vec_seq_len.begin(), vec_seq_len.end(), mean_len);
  }

  // simulate sequences
  unsigned idx_chr = 1;
  for (auto l : vec_seq_len) {
    string seq;
    generateRandomDnaSeq(seq, l, kmer_prof, rng);
    string id_chr = stringio::format("chr%d", idx_chr++);
    shared_ptr<SeqRecord> sp_rec(new SeqRecord(id_chr, "random sequence", seq));
    sp_rec->id_ref = sp_rec->id;
    this->records.push_back(sp_rec);
    // instantiate new referennce chromosome
    shared_ptr<ChromosomeReference> sp_chr(new ChromosomeReference());
    sp_chr->id = id_chr;
    sp_chr->length = seq.length();
    sp_chr->map_start_rec[0] = sp_rec;
    this->addChromosome(sp_chr);
  }
  this->num_records = num_seqs;

  return true;
}

/** Scan genome and store structural information.
  *  1) index chromosomes (start positions in genome)
  *  2) index unmasked regions (start positions in genome)
  *  3) count nucleotide frequencies
  *  4) index bp positions into buckets by nucleotide
  *  5) index bp positions into buckets by trinucleotides
  */
void GenomeReference::indexRecords() {
  // initialize data structures
  vector<unsigned> nuc_count(4, 0); // nucleotide counter
  nuc_pos = vector<vector<long> >(4); // nucleotide buckets
  unsigned cum_start = 0; // global start position
  vec_start_chr.clear(); // start positions of sequences
  boost::circular_buffer<string> trinuc(3);  // next trinucleotide to index
  // 3mer index
  map_3mer_pos.clear();
  char nucs[4] = { 'A', 'C', 'G', 'T' };
  for (auto n1 : nucs) {
    for (auto n2 : nucs) {
      for (auto n3 : nucs) {
        string key = string(1, n1) + string(1, n2) + string(1, n3);
        map_3mer_pos[key] = vector<long>();
      }
    }
  }

  vec_start_chr.push_back(cum_start);
  for (auto const & rec : records) {
    unsigned seq_len = rec->seq.length();
    // index unmasked regions (those that are not 'N')
    unsigned p = 0;
    bool is_new_region = false;
    string::const_iterator it = rec->seq.begin();
    while (it!=rec->seq.end()) {
      short nuc = nuc2idx(*it);
      // skip to next unmasked position
      while (it!=rec->seq.end() && nuc == -1) { it++; p++; nuc=nuc2idx(*it); }
      trinuc.clear(); // clear out previous trinucleotide
      // process seq until next masked position
      while (it!=rec->seq.end() && nuc > -1) {
        if (!is_new_region) {
          is_new_region = true;
          vec_start_masked.push_back(cum_start + p);
        }
        nuc_count[nuc]++;
        nuc_pos[nuc].push_back(cum_start + p);
        trinuc.push_back(string(1, *it)); // new nuc will push out oldest one from circular buffer
        if (trinuc.size() == 3) {
          string tn = trinuc[0] + trinuc[1] + trinuc[2];
          map_3mer_pos[tn].push_back(cum_start + p - 2);
        }
        it++; p++; masked_length++;
        nuc = nuc2idx(*it);
      }
      if (is_new_region) {
        vec_cumlen_masked.push_back(masked_length);
        is_new_region = false;
      }
    }
    cum_start += seq_len;
    vec_start_chr.push_back(cum_start);
  }
  length = vec_start_chr[num_records];
fprintf(stderr, "\nGenome stats:\n");
fprintf(stderr, "  records:\t\t%u\n", num_records);
fprintf(stderr, "  length:\t\t%u\n", length);
fprintf(stderr, "  length (masked):\t%u\n", masked_length);
fprintf(stderr, "Nucleotide counts:\n  A:%u\n  C:%u\n  G:%u\n  T:%u\n", nuc_count[0], nuc_count[1], nuc_count[2], nuc_count[3]);
  // calculate nucleotide frequencies (ACGT)
  double num_acgt = nuc_count[0] + nuc_count[1] + nuc_count[2] + nuc_count[3];
  nuc_freq[0] = nuc_count[0]/num_acgt;
  nuc_freq[1] = nuc_count[1]/num_acgt;
  nuc_freq[2] = nuc_count[2]/num_acgt;
  nuc_freq[3] = nuc_count[3]/num_acgt;
fprintf(stderr, "Nucleotide freqs:\n  A:%0.4f\n  C:%0.4f\n  G:%0.4f\n  T:%0.4f\n", nuc_freq[0], nuc_freq[1], nuc_freq[2], nuc_freq[3]);
}

void GenomeReference::clearIndex() {
  this->nuc_pos.clear();
  this->nuc_pos.shrink_to_fit();
  this->map_3mer_pos.clear();
}

Locus GenomeReference::getLocusByGlobalPos(long global_pos) const {
  int idx_seq = 0;
  while (global_pos >= vec_start_chr[idx_seq+1])
    idx_seq++;

  // compile locus info
  Locus *loc = new Locus();
  loc->idx_record = idx_seq;
  loc->id_ref = records[idx_seq]->id_ref;
  loc->start = global_pos - vec_start_chr[idx_seq];
  loc->end = loc->start + 1;

  return *loc;
}

Locus GenomeReference::getAbsoluteLocusMasked(double rel_pos) const {
  long double pos = rel_pos * masked_length;
  unsigned abs_pos_masked = floor(pos);
  // identify genomic region
  int idx_region = 0;
  while (abs_pos_masked > vec_cumlen_masked[idx_region])
    idx_region++;
  // calculate relative position within region
  unsigned rel_pos_in_region = abs_pos_masked;
  if (idx_region > 0)
    rel_pos_in_region -= vec_cumlen_masked[idx_region-1];
  // identify genomic sequence
  unsigned abs_start_region = vec_start_masked[idx_region];
  int idx_seq = 0;
  while (abs_start_region >= vec_start_chr[idx_seq+1])
    idx_seq++;
  // compile locus info
  Locus *loc = new Locus();
  loc->idx_record = idx_seq;
  loc->start = rel_pos_in_region + (abs_start_region - vec_start_chr[idx_seq]);
  loc->end = loc->start + 1;

  return *loc;
}

void GenomeReference::getSequence (
  const string id_chr,
  const TCoord start,
  const TCoord end,
  map<TCoord, string>& seqs
) const
{
  // perform sanity checks
  assert( chromosomes.count(id_chr) > 0 ); // chromosome exists
  assert( start < end );

  shared_ptr<ChromosomeReference> chr = chromosomes.at(id_chr);
  // do coordinates exceed chromosome limits? (should not happen)
  assert( (start >= 0) && (end <= chr->length) );

  // skip over SeqRecords located upstream of target range
  auto it_start_rec = chr->map_start_rec.begin();
  while ( (it_start_rec != chr->map_start_rec.end()) &&
          (it_start_rec->first+it_start_rec->second->seq.length() < start) )
    ++it_start_rec;

  // add SeqRecords within target range
  while ( (it_start_rec != chr->map_start_rec.end()) &&
          (it_start_rec->first < end) ) {
    ulong rec_start = it_start_rec->first;
    // determine local start and end (within current sequence)
    ulong loc_start = start-rec_start;
    ulong loc_len_max = it_start_rec->second->seq.length();
    ulong loc_len = min(end-max(rec_start, start), loc_len_max);
    seqs[loc_start] = it_start_rec->second->seq.substr(loc_start, loc_len);
    ++it_start_rec;
  }
}

bool
GenomeReference::getSequence (
  const string id_chr,
  const TCoord pos_start,
  const TCoord pos_end,
  string& seq
) const 
{
  // init return variable
  seq.clear();

  map<TCoord, string> map_start_seq;
  this->getSequence(id_chr, pos_start, pos_end, map_start_seq);
  if ( map_start_seq.size() == 0 ) // locus was not found in genome
    return false;

  // concatenate sequences for output
  for (auto start_seq : map_start_seq) {
    seq += start_seq.second;
  }

  return true;
}

} // namespace seqio