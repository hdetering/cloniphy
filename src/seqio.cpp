#include "seqio.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <stdlib.h> // for system() calls

using namespace std;

namespace seqio {

SeqRecord::SeqRecord(const string& id, const string& desc, const string& seq)
  : id(id), description(desc), seq(seq) {}

Genome::Genome(const char* filename) {
  length = 0;
  masked_length = 0;
  ploidy = 1;
  for (int i=0; i<4; i++)
    nuc_freq[i] = 0;

  // read records from file
  records = readFasta(filename);
  num_records = records.size();
  // initialize indices
  indexRecords();
}

/** Scan genome and store structural information.
  *  1) index chromosomes (start positions in genome)
  *  2) index unmasked regions (start positions in genome)
  *  3) count nucleotide frequencies
  */
void Genome::indexRecords() {
  // initialize nucleotide counter
  map<char, unsigned> nuc_count;
  nuc_count['A'] = 0;
  nuc_count['C'] = 0;
  nuc_count['G'] = 0;
  nuc_count['T'] = 0;
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
      // skip to next unmasked position
      while (it!=rec->seq.end() && *it == 'N') { it++; p++; }
      // skip to next masked position
      while (it!=rec->seq.end() && *it != 'N') {
        if (!is_new_region) {
          is_new_region = true;
          vec_start_masked.push_back(cum_start + p);
        }
        nuc_count[toupper(*it)]++;
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
fprintf(stderr, "Nucleotide counts:\n  A:%u\n  C:%u\n  G:%u\n  T:%u\n", nuc_count['A'], nuc_count['C'], nuc_count['G'], nuc_count['T']);
  // calculate nucleotide frequencies (ACGT)
  double num_acgt = nuc_count['A'] + nuc_count['C'] + nuc_count['G'] + nuc_count['T'];
  nuc_freq[0] = nuc_count['A']/num_acgt;
  nuc_freq[1] = nuc_count['C']/num_acgt;
  nuc_freq[2] = nuc_count['G']/num_acgt;
  nuc_freq[3] = nuc_count['T']/num_acgt;
fprintf(stderr, "Nucleotide freqs:\n  A:%0.4f\n  C:%0.4f\n  G:%0.4f\n  T:%0.4f\n", nuc_freq[0], nuc_freq[1], nuc_freq[2], nuc_freq[3]);
}

/** Increase the genome's ploidy to double its current value */
void Genome::duplicate() {
  ploidy *= 2;
  for (unsigned i=0; i<num_records; ++i) {
    SeqRecord orig = records[i];
    SeqRecord *dupl = new SeqRecord(orig.id, orig.description, orig.seq);
    records.push_back(*dupl);
    records[i].id += "_m";
    records[num_records+i].id += "_p";
  }
}

Locus Genome::getAbsoluteLocusMasked(double rel_pos) {
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

vector<SeqRecord> readFasta(const char *filename) {
  ifstream inputFile;
  inputFile.open(filename, ios::in);
  vector<SeqRecord> records;
  records = readFasta(inputFile);
  inputFile.close();

  return records;
}

vector<SeqRecord> readFasta(istream &input) {
  vector<SeqRecord> records;
  string header;
  string seq;
  string line;

  try {
    getline(input, header);
    while (input.good()) {
      while (stringio::safeGetline(input, line) && line[0]!='>') {
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

  return records;
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


Nuc charToNuc(const char c) {
  switch (c) {
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

char shiftNucleotide(const char base, const int offset) {
  Nuc nuc_old = charToNuc(base);
  Nuc nuc_new = static_cast<Nuc>((static_cast<int>(nuc_old) + offset) % 4); // TODO: this could be more generic (get rid of the hard-coded 4)
  return nucToChar(nuc_new);
}

} /* namespace seqio */
