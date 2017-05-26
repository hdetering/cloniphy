#include "seqio.hpp"
#include <cctype>
#include <cmath>
//#include <numeric> // std::accumulate
//#include <regex> // std::regex, std::smatch
#include <stdlib.h> // for system() calls

using namespace std;
using boost::icl::interval_map;
using boost::icl::interval;
using boost::str;

namespace seqio {

Locus::Locus() {}
Locus::Locus(string id, ulong start, ulong end)
: id_ref(id), start(start), end(end), idx_record(0) {}
Locus::~Locus() {}

SeqRecord::SeqRecord(const string id, const string desc, const string& seq)
  : id(id), description(desc), seq(seq), id_ref(id), chr_copy(0) {}
SeqRecord::~SeqRecord() {}

GenomeReference::GenomeReference() : length(0), masked_length(0) {} //, ploidy(1) {}

GenomeReference::GenomeReference(const char* filename) : length(0), masked_length(0) {
  //ploidy = 1;
  //regex rgx("$(\\w+):([0-9]+)-([0-9]+)"); // matches chromosome segment IDs
  //smatch matches;
  // initialize nucleotide freuquencies
  for (int i=0; i<4; i++)
    nuc_freq[i] = 0;

  // read records from file
  readFasta(filename, this->records);
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
    p_chr_ref->records.push_back(rec);
    // sanity check: chromosome IDs should be unique
    assert(this->chromosomes.count(p_chr_ref->id) == 0);
    this->addChromosome(p_chr_ref);
  }
  num_records = records.size();
}

void GenomeReference::addChromosome(shared_ptr<ChromosomeReference> sp_chr) {
  this->chromosomes[sp_chr->id] = sp_chr;
  this->vec_chr_id.push_back(sp_chr->id);
  this->vec_chr_len.push_back(sp_chr->length);
}

void GenomeReference::generate(
  const unsigned long total_len,
  const vector<double> nuc_freqs,
  RandomNumberGenerator<> &rng
) {
  this->generate(total_len, 1, nuc_freqs, rng);
}

/**
 * Generate random reference genome based on nucleotide frequencies.
 *
 * \param total_len total genome length including all sequences
 * \param num_chr number of chromosomes/sequences to generate
 * \param nuc_freqs nucleotide frequencies (A,C,G,T)
 * \param object to generate random numbers
 */
void GenomeReference::generate(
  const unsigned long total_len,
  const unsigned short num_chr,
  const vector<double> nuc_freqs,
  RandomNumberGenerator<> &rng
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
    string id_chr = (boost::format("chr%d") % i).str();
    shared_ptr<SeqRecord> sp_rec(new SeqRecord(id_chr, "random sequence", string(it_start, it_end)));
    this->records.push_back(sp_rec);
    // instantiate new referennce chromosome
    shared_ptr<ChromosomeReference> sp_chr(new ChromosomeReference());
    sp_chr->id = id_chr;
    sp_chr->length = it_end - it_start;
    sp_chr->records.push_back(sp_rec);
    this->addChromosome(sp_chr);
  }
  this->num_records = num_chr;
  this->length = gen_len;
}

/** generate random genome by given number of fragments, mean len, sd len, nuc freqs */
void GenomeReference::generate(
  const unsigned num_seqs,
  const unsigned long mean_len,
  const unsigned long sd_len,
  const std::vector<double> nuc_freqs,
  RandomNumberGenerator<>& rng)
{
  // set sequence lengths
  vector<unsigned long> vec_seq_len(num_seqs);
  if (sd_len > 0) { // sample from Gamma distribution
    function<double()> rlen = rng.getRandomGammaMeanSd(mean_len, sd_len);
    std::generate(vec_seq_len.begin(), vec_seq_len.end(), rlen);
    std::sort(vec_seq_len.begin(), vec_seq_len.end(), std::greater<unsigned long>());
  }
  else { // fixed length for all sequences
    std::fill(vec_seq_len.begin(), vec_seq_len.end(), mean_len);
  }

  // simulate sequences
  unsigned idx_chr = 0;
  for (auto l : vec_seq_len) {
    string seq;
    generateRandomDnaSeq(seq, l, nuc_freqs, rng);
    string id_chr = (boost::format("chr%d") % idx_chr++).str();
    shared_ptr<SeqRecord> sp_rec(new SeqRecord(id_chr, "random sequence", seq));
    sp_rec->id_ref = sp_rec->id;
    this->records.push_back(sp_rec);
    // instantiate new referennce chromosome
    shared_ptr<ChromosomeReference> sp_chr(new ChromosomeReference());
    sp_chr->id = id_chr;
    sp_chr->length = seq.length();
    sp_chr->records.push_back(sp_rec);
    this->addChromosome(sp_chr);
  }
  this->num_records = num_seqs;
}

/** Scan genome and store structural information.
  *  1) index chromosomes (start positions in genome)
  *  2) index unmasked regions (start positions in genome)
  *  3) count nucleotide frequencies
  *  4) index bp positions into buckets by nucleotide
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

void
GenomeInstance::addChromosome(
  shared_ptr<ChromosomeInstance> sp_chr,
  string id_chr)
{
  this->map_id_chr[id_chr].push_back(sp_chr);
  this->vec_chr.push_back(sp_chr);
  this->vec_chr_len.push_back(sp_chr->length);
}

void
GenomeInstance::deleteChromosome(
  shared_ptr<ChromosomeInstance> sp_chr,
  string id_chr)
{
  // remove chromosome from map
  this->map_id_chr[id_chr].erase(
    remove(
      this->map_id_chr[id_chr].begin(),
      this->map_id_chr[id_chr].end(),
      sp_chr),
    this->map_id_chr[id_chr].end());

  // rebuild indices
  this->vec_chr.clear();
  this->vec_chr_len.clear();
  for ( auto const & id_chr : this->map_id_chr ) {
    for ( auto const & sp_chr : id_chr.second ) {
      this->vec_chr.push_back(sp_chr);
      this->vec_chr_len.push_back(sp_chr->length);
    }
  }
}

void
GenomeInstance::duplicate(
  vector<seg_mod_t>& out_vec_seg_mod)
{
  vector<tuple<string, shared_ptr<ChromosomeInstance>>> vec_tpl_id_ci;
  // loop over chromosome IDs in GenomeInstance
  for (auto const & kv : this->map_id_chr) {
    // copy all chromosomes under current chromosome ID
    for (auto const & ci_old : kv.second) {
      shared_ptr<ChromosomeInstance> ci_new(new ChromosomeInstance());
      ci_new->copy(ci_old, out_vec_seg_mod);
      vec_tpl_id_ci.push_back(make_tuple(kv.first, ci_new));
    }
  }
  // add newly created ChromosomeInstances to GenomeInstance
  for (auto tpl_id_chr : vec_tpl_id_ci) {
    string id_chr;
    shared_ptr<ChromosomeInstance> sp_chr;
    tie(id_chr, sp_chr) = tpl_id_chr;
    this->addChromosome(sp_chr, id_chr);
  }
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
  unsigned abs_pos_masked = floor(rel_pos * masked_length);
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

ChromosomeReference::ChromosomeReference() : id(""), length(0) {}

ChromosomeInstance::ChromosomeInstance() : length(0) {}

ChromosomeInstance::ChromosomeInstance(ChromosomeReference ref)
//: ref_chr(shared_ptr<ChromosomeReference>(&ref)), length(ref.length) {
: length(ref.length) {
  // inititally chromosome consists of a single SegmentCopy
  SegmentCopy seg_copy(0, this->length);
  this->lst_segments.push_back(seg_copy);
}

vector<SegmentCopy> ChromosomeInstance::getSegmentCopiesAt(unsigned long ref_pos) {
  vector<SegmentCopy> res_segments;
  for (auto const & seg : this->lst_segments) {
    if ( ref_pos >= seg.ref_start && ref_pos < seg.ref_end ) {
      res_segments.push_back(seg);
    }
  }
  return res_segments;
}

vector<seg_mod_t> ChromosomeInstance::amplifyRegion(
  double start_rel,
  double len_rel,
  bool is_forward,
  bool is_telomeric
) {
  vector<seg_mod_t> vec_seg_mod;
  // iterators pointing to SegmentCopy
  list<SegmentCopy>::iterator it_seg_bkp_left, it_seg_bkp_right;
  list<SegmentCopy> lst_seg_new;
  bool is_left_bkp = false;
  // NOTE: start_rel value may be changed to match other parameters
  //len_rel = min(len_rel, is_forward ? 1.0-start_rel : start_rel);
  // perform some sanity checks
  assert( len_rel > 0.0 && len_rel <= 1.0 );
  if (is_forward)
    assert( start_rel + len_rel <= 1.0 );
  else
    assert( start_rel - len_rel >= 0.0 );
  // physical start, length of event
  unsigned long start_bp = start_rel * this->length;
  unsigned long len_bp = len_rel * this->length;
  // physical coordinates
  unsigned long bkp_start, bkp_end;
  if (is_forward) {
    bkp_start = start_bp;
    bkp_end = is_telomeric ? this->length : start_bp+len_bp;
  } else { // !is_forward
    bkp_end = start_bp;
    bkp_start = is_telomeric ? 0 : bkp_end-len_bp;
  }

  // identify first affected SegmentCopy
  unsigned long pos_bp = 0; // current position
  unsigned long seg_len = 0; // length of current SegmentCopy
  auto it_seg = this->lst_segments.begin();
  // NOTE: after break, pos_bp stores the start coordinate of the current SegmentCopy
  while (it_seg != this->lst_segments.end()) {
    seg_len = it_seg->ref_end - it_seg->ref_start;
    if (pos_bp+seg_len > bkp_start) break;
    pos_bp += seg_len;
    it_seg++;
  }
  assert( it_seg != this->lst_segments.end() );

  // breakpoint after beginning of current SegmentCopy?
  if (bkp_start >= pos_bp) {
    is_left_bkp = true;
    seg_len = it_seg->ref_end - it_seg->ref_start;
    // create new SegmentCopy
    unsigned long seg_new_start = it_seg->ref_start + (bkp_start - pos_bp);
    unsigned long seg_new_end = pos_bp+seg_len <= bkp_end ? it_seg->ref_end : seg_new_start+len_bp;
    assert( seg_new_start < seg_new_end );
    SegmentCopy seg_new(seg_new_start, seg_new_end);
    lst_seg_new.push_back(seg_new);
    vec_seg_mod.push_back(make_tuple(seg_new.id, it_seg->id, seg_new_start, seg_new_end));

    // will new SegmentCopies be inserted at first breakpoint?
    if (!is_forward ) {
      if ( bkp_start > pos_bp) { // split this SegmentCopy at breakpoint
        // head of SegmentCopy
        unsigned long head_len = bkp_start - pos_bp;
        unsigned long head_start = it_seg->ref_start;
        unsigned long head_end = it_seg->ref_start + head_len;
        SegmentCopy seg_head(head_start, head_end);
        vec_seg_mod.push_back(make_tuple(seg_head.id, it_seg->id, head_start, head_end));
        // tail of SegmentCopy
        unsigned long tail_start = head_end;
        unsigned long tail_end = it_seg->ref_end;
        SegmentCopy seg_tail(tail_start, tail_end);
        vec_seg_mod.push_back(make_tuple(seg_tail.id, it_seg->id, tail_start, tail_end));

        // replace existing SegmentCopy
        it_seg = this->lst_segments.erase(it_seg);
        this->lst_segments.insert(it_seg, seg_head);
        it_seg_bkp_left = this->lst_segments.insert(it_seg, seg_tail);
        // continue processing with new tail SegmentCopy
        it_seg = it_seg_bkp_left;
        pos_bp += head_len;
        is_left_bkp = false;
      } else { // start of SegmentCopy coincides with first breakpoint, no need to split
        // new SegmentCopies will be inserted before current one
        it_seg_bkp_left = it_seg;
      }
    }
  }

  // skip over SegmentCopies that are affected entirely
  while (it_seg != this->lst_segments.end()) {
    seg_len = it_seg->ref_end - it_seg->ref_start;
    if (pos_bp+seg_len >= bkp_end) break;
    // right breakpoint is after end of SegmentCopy
    if (!is_left_bkp) {
      SegmentCopy seg_copy(it_seg->ref_start, it_seg->ref_end);
      lst_seg_new.push_back(seg_copy);
      vec_seg_mod.push_back(make_tuple(seg_copy.id, it_seg->id, it_seg->ref_start, it_seg->ref_end));
    } else {
      is_left_bkp = false;
    }
    // proceed to next SegmentCopy
    pos_bp += seg_len;
    it_seg++;
  }

  // breakpoint before end of current SegmentCopy?
  if (it_seg != this->lst_segments.end() && bkp_end <= pos_bp+seg_len) {
    // new SegmentCopy needed only if this one does not contain left breakpoint
    // (otherwise the whole duplication has already been handled)
    if (!is_left_bkp) {
      unsigned long seg_new_start = it_seg->ref_start;
      unsigned long seg_new_end = seg_new_start + (bkp_end-pos_bp);
      assert( seg_new_start < seg_new_end );
      SegmentCopy seg_new(seg_new_start, seg_new_end);
      lst_seg_new.push_back(seg_new);
      vec_seg_mod.push_back(make_tuple(seg_new.id, it_seg->id, seg_new_start, seg_new_end));
    }
    if (is_forward ) {
      if (bkp_end < pos_bp+seg_len) { // split this SegmentCopy at breakpoint
        // head of SegmentCopy
        unsigned long head_len = bkp_end - pos_bp;
        unsigned long head_start = it_seg->ref_start;
        unsigned long head_end = it_seg->ref_start + head_len;
        SegmentCopy seg_head(head_start, head_end);
        vec_seg_mod.push_back(make_tuple(seg_head.id, it_seg->id, head_start, head_end));

        // tail of SegmentCopy
        unsigned long tail_start = head_end;
        unsigned long tail_end = it_seg->ref_end;
        SegmentCopy seg_tail(tail_start, tail_end);
        vec_seg_mod.push_back(make_tuple(seg_tail.id, it_seg->id, tail_start, tail_end));

        // replace existing SegmentCopy
        it_seg = this->lst_segments.erase(it_seg);
        this->lst_segments.insert(it_seg, seg_head);
        it_seg_bkp_right = this->lst_segments.insert(it_seg, seg_tail);
      } else { // amplification includes end of SegmentCopy
        it_seg_bkp_right = ++it_seg;
      }
    }
  }

  // insert new SegmentCopies
  if (is_forward) {
    this->lst_segments.insert(it_seg_bkp_right, lst_seg_new.begin(), lst_seg_new.end());
  } else {
    this->lst_segments.insert(it_seg_bkp_left, lst_seg_new.begin(), lst_seg_new.end());
  }

  // update ChromosomeInstance length
  this->length += len_bp;

  return vec_seg_mod;
}

vector<seg_mod_t> ChromosomeInstance::deleteRegion(
  double start_rel,
  double len_rel,
  bool is_forward,
  bool is_telomeric
) {
  vector<seg_mod_t> vec_seg_mod;
  // iterators pointing to SegmentCopoe
  list<SegmentCopy>::iterator it_rm_from, it_rm_to;
  list<SegmentCopy> lst_seg_new;
  // perform some sanity checks
  assert( len_rel > 0.0 && len_rel <= 1.0 );
  if (is_forward)
    assert( start_rel + len_rel <= 1.0 );
  else
    assert( start_rel - len_rel >= 0.0 );
  // physical start, length of event
  unsigned long start_bp = start_rel * this->length;
  unsigned long len_bp = len_rel * this->length;
  // physical coordinates
  unsigned long bkp_start, bkp_end;
  if (is_forward) {
    bkp_start = start_bp;
    bkp_end = is_telomeric ? this->length : start_bp+len_bp;
  } else { // !is_forward
    bkp_end = start_bp;
    bkp_start = is_telomeric ? 0 : bkp_end-len_bp;
  }

  // identify first affected SegmentCopy
  unsigned long pos_bp = 0; // current position
  unsigned long seg_len = 0; // length of current SegmentCopy
  auto it_seg = this->lst_segments.begin();
  // NOTE: after break, pos_bp stores the start coordinate of the current SegmentCopy
  while (it_seg != this->lst_segments.end()) {
    seg_len = it_seg->ref_end - it_seg->ref_start;
    if (pos_bp+seg_len >= bkp_start) break;
    pos_bp += seg_len;
    it_seg++;
  }
  assert( it_seg != this->lst_segments.end() );

  // current SegmentCopy will be removed
  // (if only part of SegmentCopy is affected, a new SegmentCopy will be introduced)
  it_rm_from = it_seg;

  // does deletion start after beginning of current SegmentCopy?
  if (bkp_start > pos_bp) {
    // head of SegmentCopy will be conserved
    unsigned long head_len = bkp_start - pos_bp;
    unsigned long head_start = it_seg->ref_start;
    unsigned long head_end = it_seg->ref_start + head_len;
    assert( head_start < head_end );
    SegmentCopy seg_head(head_start, head_end);
    lst_seg_new.push_back(seg_head);
    vec_seg_mod.push_back(make_tuple(seg_head.id, it_seg->id, head_start, head_end));
  }

  // skip over SegmentCopies to be removed completely
  while (it_seg != this->lst_segments.end()) {
    seg_len = it_seg->ref_end - it_seg->ref_start;
    if (pos_bp+seg_len >= bkp_end) break;
    pos_bp += seg_len;
    it_seg++;
  }

  // does deletion end before end of current SegmentCopy?
  if (bkp_end < pos_bp+seg_len) {
    // tail of SegmentCopy will be conserved
    unsigned long tail_len = pos_bp + seg_len - bkp_end;
    unsigned long tail_start = it_seg->ref_end - tail_len;
    unsigned long tail_end = it_seg->ref_end;
    assert( tail_start < tail_end );
    SegmentCopy seg_tail(tail_start, tail_end);
    lst_seg_new.push_back(seg_tail);
    vec_seg_mod.push_back(make_tuple(seg_tail.id, it_seg->id, tail_start, tail_end));
  }

  // remove affected SegmentCopies [it_rm_from, it_rm_to)
  it_rm_to = ++it_seg;
  this->lst_segments.erase(it_rm_from, it_rm_to);
  // insert new SegmentCopies
  this->lst_segments.insert(it_rm_to, lst_seg_new.begin(), lst_seg_new.end());
  // update ChromosomeInstance length
  this->length -= len_bp;

  return vec_seg_mod;
}

void ChromosomeInstance::copy(shared_ptr<ChromosomeInstance> ci_old, vector<seg_mod_t>& out_vec_seg_mod) {
  this->length = ci_old->length;
  for (auto const & seg_old : ci_old->lst_segments) {
    SegmentCopy seg_new(seg_old.ref_start, seg_old.ref_end);
    this->lst_segments.push_back(seg_new);
    out_vec_seg_mod.push_back(make_tuple(seg_new.id, seg_old.id, seg_old.ref_start, seg_old.ref_end));
  }
}

ostream& operator<<(ostream& lhs, const ChromosomeInstance& ci) {
  lhs << "    ChromosomeInstance<length=" << ci.length << ">" << endl;
  for (auto const & seg : ci.lst_segments) {
    lhs << "      " << seg;
  }
  return lhs;
}

GenomeInstance::GenomeInstance() {}

GenomeInstance::GenomeInstance(GenomeReference g_ref) {
  for (auto const & kv : g_ref.chromosomes) {
    ChromosomeReference chr_ref = *(kv.second);
    // initial genome state is diploid -> generate two instances of each chromosome
    shared_ptr<ChromosomeInstance> sp_chr_inst1(new ChromosomeInstance(chr_ref));
    this->vec_chr.push_back(sp_chr_inst1);
    this->vec_chr_len.push_back(sp_chr_inst1->length);
    shared_ptr<ChromosomeInstance> sp_chr_inst2(new ChromosomeInstance(chr_ref));
    this->vec_chr.push_back(sp_chr_inst2);
    this->vec_chr_len.push_back(sp_chr_inst2->length);
    // sanity check: chromsome IDs should be unique
    assert(this->map_id_chr.count(chr_ref.id)==0);
    // shared_ptr<ChromosomeInstance> sp_chr_inst1(up_chr_inst1);
    // shared_ptr<ChromosomeInstance> sp_chr_inst2(up_chr_inst2);
    this->map_id_chr[chr_ref.id] = { sp_chr_inst1, sp_chr_inst2 };
  }
}

vector<SegmentCopy>
GenomeInstance::getSegmentCopiesAt (
  string id_chr,
  unsigned long ref_pos)
{
  vector<SegmentCopy> res_segments;
  // perform sanity checks
  assert( this->map_id_chr.count(id_chr) > 0 );

  // get SegmentCopies from all ChromosomeInstances
  for (auto const sp_chr : this->map_id_chr[id_chr]) {
    vector<SegmentCopy> chr_segs = sp_chr->getSegmentCopiesAt(ref_pos);
    res_segments.insert(res_segments.end(), chr_segs.begin(), chr_segs.end());
  }

  return res_segments;
}

void
GenomeInstance::getCopyNumberStates(
  map<unsigned, vector<shared_ptr<Locus>>>& map_cn_loci) {
  // for all chromosomes
  for ( auto const & id_chr : this->map_id_chr ) {
    string id = id_chr.first;
    interval_map<unsigned long, int> imap_reg_cn;
    // for all ChromosomeInstances
    for ( auto const & sp_chr : id_chr.second ) {
      // for all SegmentCopies
      for ( auto const & seg : sp_chr->lst_segments ) {
        auto i_reg = interval<unsigned long>::right_open(seg.ref_start, seg.ref_end);
        imap_reg_cn += make_pair(i_reg, 1);
      }
    }

    // collect genomic regions and CN states
    for (auto const & i_reg_cn : imap_reg_cn) {
      auto reg = i_reg_cn.first;
      int cn = i_reg_cn.second;
      if (map_cn_loci.count(cn) == 0)
        map_cn_loci[cn] = vector<shared_ptr<Locus>>();
      shared_ptr<Locus> sp_loc(new Locus(id, reg.lower(), reg.upper()));
      map_cn_loci[cn].push_back(sp_loc);
    }
  }
}

ostream& operator<<(ostream& lhs, const GenomeInstance& gi) {
  lhs << "GenomeInstance" << endl;
  lhs << "--------------" << endl;
  for (auto const & id_chr : gi.map_id_chr) {
    unsigned num_ci = id_chr.second.size();
    lhs << "  Chromosome \"" << id_chr.first << "\" (" << num_ci << " instances):" << endl;
    for (auto const & ci : id_chr.second) {
      lhs << *ci;
    }
  }
  return lhs;
}

SegmentCopy::SegmentCopy()
: id(boost::uuids::random_generator()()) {}

SegmentCopy::~SegmentCopy() {}

SegmentCopy::SegmentCopy(unsigned long start, unsigned long end)
: id(boost::uuids::random_generator()()), ref_start(start), ref_end(end) {}

ostream& operator<<(ostream& lhs, const SegmentCopy& seg) {
  lhs << "SegmentCopy<uuid=" << seg.id << "> ";
  lhs << "[" << seg.ref_start << ", " << seg.ref_end << ")";
  lhs << endl;
  return lhs;
}


/*------------------------------------*/
/*           Utility methods          */
/*------------------------------------*/

void readFasta(const char *filename, vector<shared_ptr<SeqRecord>> &records) {
  ifstream inputFile;
  inputFile.open(filename, ios::in);
  readFasta(inputFile, records);
  inputFile.close();
}

void readFasta(istream &input, vector<shared_ptr<SeqRecord>> &records) {
  string header;
  string seq;
  string line;

  try {
    stringio::safeGetline(input, header);
    while (input.good()) {
      while (stringio::safeGetline(input, line)) {
        if (line.length()>0 && line[0]=='>')
          break;
        seq += line;
      }
      size_t space_pos = header.find(' ');
      string seq_id = header.substr(1, space_pos-1);
      string seq_desc = header.substr(space_pos+1);
      //SeqRecord rec = {seq_id, seq_desc, seq};
      shared_ptr<SeqRecord> sp_rec(new SeqRecord(seq_id, seq_desc, seq));
      // get properties from description
      vector<string> desc_parts = stringio::split(seq_desc, ';');
      for (string part : desc_parts) {
        vector<string> kv = stringio::split(part, '=');
        if (kv.size() == 2)
          if (kv[0] == "id_ref")
            sp_rec->id_ref = kv[1];
      }
      records.push_back(sp_rec);
      header = line;
      seq = "";
    }
  } catch(exception e) {
    cerr << e.what() << endl;
    cerr << "Possibly premature end of FASTA file?" << endl;
  }
}

/** Write SeqRecords to FASTA file, using a defined line width. */
int writeFasta(
  const vector<shared_ptr<SeqRecord>>& sequences,
  const string filename,
  int line_width)
{
  int num_records = 0;
  ofstream ofs;
  ofs.open(filename);
  num_records = writeFasta(sequences, ofs, line_width);
  ofs.close();

  return num_records;
}

/** Write SeqRecords to FASTA file, using a defined line width. */
int writeFasta(
  const vector<shared_ptr<SeqRecord>>& seqs,
  ostream &output,
  int line_width)
{
  int recCount = 0;
  for (auto const & rec : seqs) {
    // seq description confuses ART (mismatch of SAM header with REF field)
    //output << str(boost::format(">%s id_ref=%s\n") % rec.id % rec.id_ref);
    output << str(boost::format(">%s\n") % rec->id);
    string::const_iterator it_seq = rec->seq.begin();
    while (it_seq != rec->seq.end()) {
      for (int i=0; i<line_width && it_seq!=rec->seq.end(); ++i)
        output << *it_seq++;
      output << endl;
    }
    recCount++;
  }
  return recCount;
}

void indexFasta(const char *filename) {
  // TODO: implement this function
  string fn(filename);
  string cmd = "samtools faidx " + fn;
  system(cmd.c_str());
}

unsigned long generateRandomDnaSeq(
  string &seq,
  const unsigned long total_len,
  const vector<double>nuc_freqs,
  RandomNumberGenerator<> &rng
) {
  function<short()> ridx = rng.getRandomIndexWeighted(nuc_freqs);

  // generate random sequence
  unsigned long gen_len = 0;
  while (gen_len < total_len) {
    // pick random nucleotide
    char nuc = idx2nuc(ridx());
    // add to genomic sequence
    seq += nuc;
    gen_len++;
  }
  return gen_len;
}

/** Simulate allelic dropout (ADO) events
  *   percentage: fraction of genome to mask with 'N's
  *   frag_len:   length of masked fragments (runs of 'N's) */
void simulateADO_old(const string fn_input,
                 const float percentage,
                 const int frag_len,
                 function<double()>& random) {
fprintf(stderr, "Simulating ADO for file '%s'\n", fn_input.c_str());
  // read input genome
  GenomeReference g = GenomeReference(fn_input.c_str());

  // calculate number of masked regions needed
  int num_frags = (int)(ceil((g.masked_length * percentage) / frag_len));
  // pick locations for masked regions
  // (multiples of fragment length so they don't overlap)
  // TODO: avoid a start position to be picked twice
  vector<double> vec_pos_frag;
  for (int i=0; i<num_frags; ++i) {
    float r = random();
    double pos_folded = (unsigned)trunc(r*(g.masked_length/frag_len));
    double pos_abs = pos_folded * frag_len;
    vec_pos_frag.push_back(pos_abs/g.masked_length);
  }
  sort(vec_pos_frag.begin(), vec_pos_frag.end());

  // perform actual masking
  for (size_t i=0; i<vec_pos_frag.size(); ++i) {
    Locus loc = g.getAbsoluteLocusMasked(vec_pos_frag[i]);
//fprintf(stderr, "\tmasking: %s:%u-%u\n", g.records[loc.idx_record].id.c_str(), loc.start, loc.start+1000);
    for (int p=0; p<frag_len; ++p)
      g.records[loc.idx_record]->seq[loc.start+p] = 'N';
  }

  // create output filename
  size_t pos_dot = fn_input.find_last_of(".");
  string fn_output = fn_input;
  fn_output.insert(pos_dot, ".ado");

  // write modified genome to file
fprintf(stderr, "Writing modified genome to file '%s'\n", fn_output.c_str());
  ofstream f_out;
  f_out.open(fn_output.c_str());
  writeFasta(g.records, f_out);
  f_out.close();
}

/** Simulate allelic dropout (ADO) events
  *   percentage: fraction of genome to mask with 'N's
  *   frag_len:   length of masked fragments (runs of 'N's) */
void simulateADO(const string fn_input,
                 const unsigned genome_len,
                 const float percentage,
                 const int frag_len,
                 function<double()>& random) {
fprintf(stderr, "Simulating ADO for file '%s'\n", fn_input.c_str());
  // TODO: might be worth doing this in streaming fashion
  // open genome input file
  //ifstream f_in;
  //f_in.open(fn_input);
  GenomeReference g = GenomeReference(fn_input.c_str());

  // create output filenames
  size_t pos_dot = fn_input.find_last_of(".");
  string fn_output = fn_input;
  string fn_bed = fn_input;
  fn_output.insert(pos_dot, ".ado");
  fn_bed.replace(pos_dot, string::npos, ".ado.bed");

  // generate BED file keeping track of ADO fragments
  ofstream f_bed;
  f_bed.open(fn_bed.c_str());

  // calculate number of masked regions needed
  int num_frags = (int)(ceil((genome_len * percentage) / frag_len));
  double r_frag = ((double)num_frags) / genome_len; // fragment start rate
  // pick locations for masked regions
  // (multiples of fragment length so they don't overlap)
  // TODO: avoid a start position to be picked twice
  bool is_fragment = false;
  int idx_char = 0;
  for (unsigned idx_chr=0; idx_chr<g.records.size(); ++idx_chr) {
    unsigned seq_len = g.records[idx_chr]->seq.length();
    for (unsigned idx_nuc=0; idx_nuc<seq_len; ++idx_nuc) {
      if (!is_fragment) {
        float r = random();
        if (r <= r_frag) {
          is_fragment = true;
          f_bed << boost::format("%s\t%u\t%u\n") % g.records[idx_chr]->id % idx_nuc % min(idx_nuc+frag_len, seq_len);
        }
      }
      if (is_fragment) {
        if (idx_char < frag_len) {
          g.records[idx_chr]->seq[idx_nuc] = 'N';
          idx_char++;
        }
        else {
          is_fragment = false;
          idx_char = 0;
        }
      }
    }
  }
  f_bed.close();

  // write modified genome to file
fprintf(stderr, "Writing modified genome to file '%s'\n", fn_output.c_str());
  ofstream f_out;
  f_out.open(fn_output.c_str());
  writeFasta(g.records, f_out);
  f_out.close();
}


Nuc charToNuc(const char c) {
  switch (toupper(c)) {
    case 'A':
      return A;
    case 'C':
      return C;
    case 'G':
      return G;
    case 'T':
      return T;
    default:
      return N;
  }
}

char nucToChar(const Nuc n) {
  switch (n) {
    case A:
      return 'A';
    case C:
      return 'C';
    case G:
      return 'G';
    case T:
      return 'T';
    default:
      return 'N';
  }
}

string nucToString(const Nuc n) {
  stringstream ss;
  string s;

  ss << nucToChar(n);
  ss >> s;

  return s;
}

// TODO: deprecated?
char shiftNucleotide(const char base, const int offset) {
  Nuc nuc_old = charToNuc(base);
  Nuc nuc_new = static_cast<Nuc>((static_cast<int>(nuc_old) + offset) % 4); // TODO: this could be more generic (get rid of the hard-coded 4)
  return nucToChar(nuc_new);
}

} /* namespace seqio */
