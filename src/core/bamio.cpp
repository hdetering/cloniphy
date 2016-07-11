#include "bamio.hpp"
using boost::format;
using boost::str;
using seqan::BamAlignmentRecord;
using seqan::BamFileIn;
using seqan::BamFileOut;
using seqan::BamHeader;
using seqan::CharString;
using seqan::FormattedFileContext;
using seqan::getAbsolutePath;
using namespace std;

namespace bamio {

void mutateReads(
  string fn_out,
  string fn_in,
  int num_mutations,
  vector<vario::Variant> variants,
  treeio::Tree<Clone> tree,
  RandomNumberGenerator<> rng)
{
  // extract info about visible clones
  vector<Clone*> vec_vis_clones = tree.getVisibleNodes();
  vector<int> vec_clone_idx;
  vector<string> vec_clone_lbl;
  vector<double> vec_clone_weight;
  for (Clone *c : vec_vis_clones) {
    fprintf(stdout, "clone %d: \"%s\"\t(%.4f)\n", c->index, c->label.c_str(), c->weight);
    vec_clone_idx.push_back(c->index);
    vec_clone_lbl.push_back(c->label);
    vec_clone_weight.push_back(c->weight);
  }

  // get mutation matrix
  int num_nodes = tree.m_numNodes;
  vector<vector<bool>> mm(num_nodes, vector<bool>(num_mutations, false));
  tree.m_root->populateMutationMatrixRec(mm);

  // postprocessing: change CHR field to match diploid identifiers (used in BAM file)
  vector<Variant> var_sorted = Variant::sortByPositionPoly(variants);
  for (Variant &var : var_sorted)
    var.chr += str(format("_%d") % var.chr_copy);

  /*------------------*
   * process BAM file *
   *------------------*/

  function<int()> rand_idx = rng.getRandomIndexWeighted(vec_clone_weight);
  CharString bamFileName = getAbsolutePath(fn_in.c_str());
  BamFileIn bamFileIn;
  if (!open(bamFileIn, toCString(bamFileName))) {
    std::cerr << "ERROR: Could not open " << bamFileName << std::endl;
    return;
  };

  // the BAM file context gives access to reference names etc.
  //typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;
  typedef seqan::StringSet<CharString> TNameStore;
  typedef seqan::NameStoreCache<TNameStore> TNameStoreCache;
  typedef seqan::BamIOContext<TNameStore, TNameStoreCache, seqan::Dependent<> > TBamContext;
  TBamContext & bamContextIn = context(bamFileIn);

  // create output SAM file for mutated reads
  ofstream fs_bamout;
  fs_bamout.open("bulk_reads.mutated.sam");
  BamFileOut bamFileOut(context(bamFileIn), fs_bamout, seqan::Sam());
  //BamFileOut bamFileOut(context(bamFileIn), cout, seqan::Sam());
  // NOTE: reading the header first is MANDATORY!
  BamHeader headerIn, headerOut;
  readHeader(headerIn, bamFileIn);
  TNameStore refNames; // keep haploid references
  TNameStore refNamesMap; // map diploid references to haploid names
  typedef seqan::Member<TBamContext, seqan::LengthStoreMember>::Type TLengthStore;
  TLengthStore refLengths;
  seqan::BamHeaderRecord headRec;
  for (auto headRec : headerIn) {
    // include only haploid references (seq sim was done on diploid genome)
    if (headRec.type == seqan::BAM_HEADER_REFERENCE) {
      CharString seqName, seqLen;
      getTagValue(seqName, "SN", headRec);
      getTagValue(seqLen, "LN", headRec);
      string seq_name(toCString(seqName));
      int pos_delim = seq_name.rfind('_');
      if (pos_delim == string::npos) {
        fprintf(stderr, "[WARN] invalid sequence name: %s (missing copy info).\n", seq_name.c_str());
        appendValue(refNames, seqName);
        appendValue(refLengths, atoi(toCString(seqLen)));
        appendValue(headerOut, headRec);
      }
      else {
        CharString seqNameHap(seq_name.substr(0, pos_delim));
        //appendValue(refNames, seqNameHap);
        appendValue(refNamesMap, seqNameHap);
        int copy = atoi(seq_name.substr(pos_delim+1).c_str());
        if (copy == 0) {
          appendValue(refNames, seqNameHap);
          appendValue(refLengths, atoi(toCString(seqLen)));
          setTagValue("SN", seqNameHap, headRec);
          appendValue(headerOut, headRec);
        }
      }
    }
    else {
      appendValue(headerOut, headRec);
    }
  }
  // Context to map diploid reads to haploid references
  TNameStoreCache refNamesMapCache(refNamesMap);
  TBamContext bamContextMap(refNamesMap, refNamesMapCache);
  // Context (haploid) for output BAM header
  TNameStoreCache refNamesCache(refNames);
  TBamContext bamContextOut(refNames, refNamesCache);
  setContigLengths(bamContextOut, refLengths);

  //writeHeader(bamFileOut, headerOut);
  write(bamFileOut.iter, headerOut, bamContextOut, bamFileOut.format);

  ofstream fs_fq, fs_log;
  fs_fq.open(fn_out);
  fs_log.open("bamio_bulk.log");
  BamAlignmentRecord read1, read2;
  while (!atEnd(bamFileIn)) {
    readRecord(read1, bamFileIn);
    readRecord(read2, bamFileIn);
    // pick random clone (weighted) to which next read pair belongs
    int c_idx = vec_clone_idx[rand_idx()];
    int r1_begin = read1.beginPos;
    int r2_begin = read2.beginPos;
    int r1_len = getAlignmentLengthInRef(read1);
    int r2_len = getAlignmentLengthInRef(read2);
    CharString r1_ref = contigNames(bamContextIn)[read1.rID];
    CharString r2_ref = contigNames(bamContextIn)[read2.rID];
    char r1_rc = seqan::hasFlagRC(read1) ? '+' : '-';
    char r2_rc = seqan::hasFlagRC(read2) ? '+' : '-';
    string r1_qual = toCString(read1.qual);
    string r2_qual = toCString(read2.qual);
    fs_log << str(format("%s:%s(%d,%d) -- %s:%s(%d,%d)\n")
      % r1_ref % r1_rc % r1_begin % (r1_begin+r1_len)
      % r2_ref % r2_rc % r2_begin % (r2_begin+r2_len));

    // modify read pair to match assigned clone
    read1.qName += str(format("-%d/1") % c_idx);
    read2.qName += str(format("-%d/2") % c_idx);

    // mutate reads

    auto it_var = var_sorted.begin();
    // advance variant iterator to position in current chromosome
    while (it_var!=var_sorted.end() && it_var->chr!=toCString(r1_ref))
      ++it_var;
    // advance variant iterator to first position past begin of first read
    while (it_var!=var_sorted.end() && it_var->pos<r1_begin && (it_var+1)->chr==toCString(r1_ref))
      ++it_var;
    // identify variables affecting read pair
    while (it_var!=var_sorted.end() && it_var->pos<r2_begin+r2_len && it_var->chr==toCString(r1_ref)) {
      if (!mm[c_idx][it_var->idx_mutation]) { // does read pair's clone carry mutation?
        ++it_var;
        continue;
      }
      int r1_var_pos = it_var->pos - r1_begin;
      if (r1_var_pos >= 0 && r1_var_pos < r1_len) { // mutate read1
        fs_log << str(format("%s:%d\t%s->%s\n") % toCString(read1.qName) % r1_var_pos % it_var->alleles[0].c_str() % it_var->alleles[1].c_str());
        read1.seq[r1_var_pos] = it_var->alleles[1][0];
      }
      int r2_var_pos = it_var->pos - r2_begin;
      if (r2_var_pos >= 0 && r2_var_pos < r2_len) { // mutate read2
        fs_log << str(format("%s:%d\t%s->%s\n") % toCString(read2.qName) % r2_var_pos % it_var->alleles[0].c_str() % it_var->alleles[1].c_str());
        read2.seq[r2_var_pos] = it_var->alleles[1][0];
      }
      ++it_var;
    }

    // BAM output

    //CharString r1_ref_hap = refNames[read1.rID];
    //CharString r2_ref_hap = refNames[read2.rID];
    write(bamFileOut.iter, read1, bamContextMap, bamFileOut.format);
    write(bamFileOut.iter, read2, bamContextMap, bamFileOut.format);

    // FASTQ output

    // re-reverse complement reads on opposite strand (to avoid strand bias)
    if (seqan::hasFlagRC(read1)) {
      seqan::reverseComplement(read1.seq);
      reverse(r1_qual.begin(), r1_qual.end());
    }
    if (seqan::hasFlagRC(read2)) {
      seqan::reverseComplement(read2.seq);
      reverse(r2_qual.begin(), r2_qual.end());
    }

    // write (potentially modified) reads to FASTQ file
    fs_fq << str(format("@%s\n%s\n+\n%s\n") % read1.qName % read1.seq % r1_qual);
    fs_fq << str(format("@%s\n%s\n+\n%s\n") % read2.qName % read2.seq % r2_qual);

    //fs_idx << r1_begin << " (" << r1_len << ")" <<  "  " << r2_begin << " (" << r2_len << ")" << endl;
    //fprintf(stderr, "%d\n", c_idx);
  }

  fs_bamout.close();
  fs_fq.close();
  fs_log.close();

}


} // namespace bamio
