#include "seqio.hpp"
#include <cctype>
#include <cmath>

using namespace std;

namespace seqio {

/*------------------------------------*/
/*           Utility methods          */
/*------------------------------------*/

/* DEPRECATED
void readFasta (
  vector<shared_ptr<SeqRecord>> &records, 
  const char *filename
) {
  ifstream inputFile;
  inputFile.open(filename, ios::in);
  readFasta(records, inputFile);
  inputFile.close();
} */

void readFasta (
  vector<shared_ptr<SeqRecord>> &records, 
  const char *filename,
  const bool use_whitelist,
  const TLocusMap& map_seq_loci
) {
  ifstream inputFile;
  inputFile.open(filename, ios::in);
  readFasta(records, inputFile, use_whitelist, map_seq_loci);
  inputFile.close();
}

void readFasta (
  vector<shared_ptr<SeqRecord>> &records, 
  istream &input,
  const bool use_whitelist,
  const TLocusMap& map_loci 
) {
  string header;
  string seq;
  string line;

  try {
    stringio::safeGetline(input, header);
    while (input.good()) {
      // read sequence data
      while (stringio::safeGetline(input, line)) {
        if (line.length()>0 && line[0]=='>')
          break;
        seq += line;
      }
      
      // parse header line
      size_t space_pos = header.find(' ');
      string seq_id = header.substr(1, space_pos-1);
      string seq_desc = header.substr(space_pos+1);
      
      // check if record in list of sequences to import
      if (map_loci.find(seq_id) != map_loci.end() || !use_whitelist) {
        shared_ptr<SeqRecord> sp_rec(new SeqRecord(seq_id, seq_desc, seq));
        // get properties from description
        vector<string> desc_parts = stringio::split(seq_desc, ';');
        for (string part : desc_parts) {
          vector<string> kv = stringio::split(part, '=');
          if (kv.size() == 2)
            if (kv[0] == "id_ref")
              sp_rec->id_ref = kv[1];
        }
        // add record to output
        records.push_back(sp_rec);
      }

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
    output << stringio::format(">%s\n", rec->id.c_str()).c_str();
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

unsigned long 
generateRandomDnaSeq (
  string &seq,
  const unsigned long total_len,
  const vector<double> nuc_freqs,
  RandomNumberGenerator &rng
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

unsigned long 
generateRandomDnaSeq (
  string &seq,
  const unsigned long total_len,
  const KmerProfile& kmer_prof,
  RandomNumberGenerator &rng
) {
  unsigned idx;
  string kmer;

  // pick first kmer at random
  idx = rng.getRandomIndexWeighted(kmer_prof.m_vec_weight)();
  kmer = kmer_prof.m_vec_kmers[idx];
  seq = kmer;
  unsigned long gen_len = kmer.length();

  // generate random sequence
  while (gen_len < total_len) {
    // get last k-1 characters from previous kmer
    string sfx = kmer.substr(1);
    // pick random kmer starting with last k-1 characters
    idx = rng.getRandomIndexWeighted(kmer_prof.m_idx_pfx.at(sfx).m_vec_weight)();
    kmer = kmer_prof.m_idx_pfx.at(sfx).m_vec_kmers[idx];
    // add last character in kmer to genomic sequence
    seq += kmer.substr(kmer.length()-1);
    gen_len++;
  }
  return gen_len;
}

/** Simulate allelic dropout (ADO) events
  *   percentage: fraction of genome to mask with 'N's
  *   frag_len:   length of masked fragments (runs of 'N's) */
/* void simulateADO_old(const string fn_input,
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
 */
/** Simulate allelic dropout (ADO) events
  *   percentage: fraction of genome to mask with 'N's
  *   frag_len:   length of masked fragments (runs of 'N's) */
/* void simulateADO(const string fn_input,
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
          f_bed << stringio::format("%s\t%u\t%u\n", (g.records[idx_chr]->id).c_str(), 
                                    idx_nuc, min(idx_nuc+frag_len, seq_len)).c_str();
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
} */


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
