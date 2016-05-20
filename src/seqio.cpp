#include "seqio.hpp"
#include <algorithm>
#include <boost/format.hpp>
#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdlib.h> // for system() calls

using namespace std;

namespace seqio {

SeqRecord::SeqRecord(const string id, const string desc, const string& seq)
  : id(id), description(desc), seq(seq), id_ref(id), copy(0) {}

Genome::Genome(const char* filename) {
  length = 0;
  masked_length = 0;
  ploidy = 1;
  for (int i=0; i<4; i++)
    nuc_freq[i] = 0;

  // read records from file
  readFasta(filename, records);
  num_records = records.size();
  // initialize indices
  indexRecords();
}

/** Scan genome and store structural information.
  *  1) index chromosomes (start positions in genome)
  *  2) index unmasked regions (start positions in genome)
  *  3) count nucleotide frequencies
  *  4) index bp positions into buckets by nucleotide
  */
void Genome::indexRecords() {
  // initialize nucleotide counter
  vector<unsigned> nuc_count(4, 0);
  // initialize nucleotide buckets
  nuc_pos = vector<vector<long> >(4);
  // initialize global start position
  unsigned cum_start = 0;
  vec_start_chr.push_back(cum_start);
  for (vector<SeqRecord>::const_iterator rec=records.begin(); rec!=records.end(); ++rec) {
    unsigned seq_len = rec->seq.length();
    // index unmasked regions (those that are not 'N')
    unsigned p = 0;
    bool is_new_region = false;
    string::const_iterator it = rec->seq.begin();
    while (it!=rec->seq.end()) {
      short nuc = nuc2idx(*it);
      // skip to next unmasked position
      while (it!=rec->seq.end() && nuc > -1) { it++; p++; }
      // skip to next masked position
      while (it!=rec->seq.end() && nuc > -1) {
        if (!is_new_region) {
          is_new_region = true;
          vec_start_masked.push_back(cum_start + p);
        }
        //nuc_count[toupper(*it)]++;

        nuc_count[nuc]++;
        //map_nuc_pos[toupper(*it)].push_back(cum_start + p);
        nuc_pos[nuc].push_back(cum_start + p);
        it++; p++; masked_length++;
      }
      if (is_new_region)
        vec_cumlen_masked.push_back(masked_length);
        is_new_region = false;
    }
    cum_start += seq_len;
    vec_start_chr.push_back(cum_start);
  }
  length = vec_start_chr[num_records];
fprintf(stderr, "\nGenome stats:\n\trecords:\t%u\n\tlength:\t%u\n", num_records, length);
fprintf(stderr, "Nucleotide counts:\n  A:%u\n  C:%u\n  G:%u\n  T:%u\n", nuc_count[0], nuc_count[1], nuc_count[2], nuc_count[3]);
  // calculate nucleotide frequencies (ACGT)
  double num_acgt = nuc_count[0] + nuc_count[1] + nuc_count[2] + nuc_count[3];
  nuc_freq[0] = nuc_count[0]/num_acgt;
  nuc_freq[1] = nuc_count[1]/num_acgt;
  nuc_freq[2] = nuc_count[2]/num_acgt;
  nuc_freq[3] = nuc_count[3]/num_acgt;
fprintf(stderr, "Nucleotide freqs:\n  A:%0.4f\n  C:%0.4f\n  G:%0.4f\n  T:%0.4f\n", nuc_freq[0], nuc_freq[1], nuc_freq[2], nuc_freq[3]);
}

/** Increase the genome's ploidy to double its current value */
void Genome::duplicate() {
  ploidy *= 2;
  for (unsigned i=0; i<num_records; ++i) {
    SeqRecord orig = records[i];
    SeqRecord *dupl = new SeqRecord(orig.id, orig.description, orig.seq);
    dupl->copy = orig.copy+1;
    records.push_back(*dupl);
    records[i].id += "_m";
    records[num_records+i].id += "_p";
  }
}

Locus Genome::getAbsoluteLocusMasked(double rel_pos) const {
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
  loc->length = 1;

  return *loc;
}

/*------------------------------------*/
/*           Utility methods          */
/*------------------------------------*/

void readFasta(const char *filename, vector<SeqRecord> &records) {
  ifstream inputFile;
  inputFile.open(filename, ios::in);
  readFasta(inputFile, records);
  inputFile.close();
}

void readFasta(istream &input, vector<SeqRecord> &records) {
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
      SeqRecord rec(seq_id, seq_desc, seq);
      records.push_back(rec);
      header = line;
      seq = "";
    }
  } catch(exception e) {
    cerr << e.what() << endl;
    cerr << "Possibly premature end of FASTA file?" << endl;
  }
}

/** Write SeqRecords to FASTA file, using a defined line width. */
int writeFasta(const vector<SeqRecord> &seqs, ostream &output, int line_width) {
  int recCount = 0;
  for (vector<SeqRecord>::const_iterator rec=seqs.begin(); rec!=seqs.end(); ++rec) {
    output << ">" << rec->id << endl;
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

/** Simulate allelic dropout (ADO) events
  *   percentage: fraction of genome to mask with 'N's
  *   frag_len:   length of masked fragments (runs of 'N's) */
void simulateADO_old(const string fn_input,
                 const float percentage,
                 const int frag_len,
                 boost::function<double()>& random) {
fprintf(stderr, "Simulating ADO for file '%s'\n", fn_input.c_str());
  // read input genome
  Genome g = Genome(fn_input.c_str());

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
      g.records[loc.idx_record].seq[loc.start+p] = 'N';
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
                 boost::function<double()>& random) {
fprintf(stderr, "Simulating ADO for file '%s'\n", fn_input.c_str());
  // TODO: might be worth doing this in streaming fashion
  // open genome input file
  //ifstream f_in;
  //f_in.open(fn_input);
  Genome g = Genome(fn_input.c_str());

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
    unsigned seq_len = g.records[idx_chr].seq.length();
    for (unsigned idx_nuc=0; idx_nuc<seq_len; ++idx_nuc) {
      if (!is_fragment) {
        float r = random();
        if (r <= r_frag) {
          is_fragment = true;
          f_bed << boost::format("%s\t%u\t%u\n") % g.records[idx_chr].id % idx_nuc % min(idx_nuc+frag_len, seq_len);
        }
      }
      if (is_fragment) {
        if (idx_char < frag_len) {
          g.records[idx_chr].seq[idx_nuc] = 'N';
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
