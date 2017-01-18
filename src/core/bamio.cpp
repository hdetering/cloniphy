#include "bamio.hpp"
#include <cstdio> // std::remove
#include <cstdlib> // system()
#include <map>
using boost::str;
using namespace std;
using namespace seqan;

namespace bamio {

ArtWrapper::ArtWrapper(string path) : bin_path(path) {
  // set default values
  read_len = 100;
  frag_len_mean = 500;
  frag_len_sd = 20;
  fold_cvg = 50;
  out_pfx = "art_sim";
  out_sam = true;
  out_aln = false;
  seq_sys = "HS25";
  fn_ref_fa = "";
  do_keep_fq = false;
}

int ArtWrapper::run(string out_pfx) {
  this->out_pfx = out_pfx;

  // build ART command line
  string art_cmd = bin_path;
  art_cmd += str(boost::format(" -l %d") % read_len);
  art_cmd += str(boost::format(" -p -m %d -s %d") % frag_len_mean % frag_len_sd);
  art_cmd += str(boost::format(" -f %.2f") % fold_cvg);
  art_cmd += " -ss " + seq_sys;
  art_cmd += " -i " + fn_ref_fa;
  art_cmd += " -o " + out_pfx;
  art_cmd += (out_sam ? " -sam" : "");
  art_cmd += (out_aln ? "" : " -na");

  fprintf(stderr, "---\nRunning ART with the following paramters:\n%s\n", art_cmd.c_str());
  int res_art = system(art_cmd.c_str());
  if (res_art != EXIT_SUCCESS) {
    fprintf(stderr, "[ERROR] ART call had non-zero return value.\n");
    return EXIT_FAILURE;
  }

  if (!do_keep_fq) {
    remove((out_pfx + "1.fq").c_str());
    remove((out_pfx + "2.fq").c_str());
  }

  return res_art;
}

/** Spike in germline mutations to SAM input. */
void mutateReads(
  const string fn_sam_out,
  const string fn_sam_in,
  const vario::VariantSet &variants,
  const short ploidy,
  RandomNumberGenerator<> &rng)
{
  // make sure input file exists
  //CharString bamFileName = getAbsolutePath(fn_sam_in.c_str());
  string bamFileName = fn_sam_in;
  BamFileIn bamFileIn;
  if (!open(bamFileIn, toCString(bamFileName))) {
    std::cerr << "ERROR: Could not open " << bamFileName << std::endl;
    return;
  };

  // the BAM file context gives access to reference names etc.
  typedef seqan::StringSet<CharString> TNameStore;
  typedef seqan::NameStoreCache<TNameStore> TNameStoreCache;
  typedef seqan::BamIOContext<TNameStore, TNameStoreCache, seqan::Dependent<> > TBamContext;
  TBamContext & bamContextIn = context(bamFileIn);

  // create output SAM file for mutated reads
  ofstream fs_bamout;
  fs_bamout.open(fn_sam_out.c_str());
  BamFileOut bamFileOut(context(bamFileIn), fs_bamout, seqan::Sam());
  // NOTE: reading the header first is MANDATORY! (advances file pointer)
  BamHeader header;
  readHeader(header, bamFileIn);
  writeHeader(bamFileOut, header);

  // initialize random function
  function<short()> rphase = rng.getRandomFunctionInt(short(0), short(ploidy-1));

  // store coverage depth and variant allele frequency for each variant
  map<string, unsigned> map_var_cvg;
  map<string, unsigned> map_var_vaf;
  for (Variant var : variants.vec_variants) {
    map_var_cvg[var.id] = 0;
    map_var_vaf[var.id] = 0;
  }

  // process reads
  BamAlignmentRecord read1, read2;
  while (!atEnd(bamFileIn)) {
    readRecord(read1, bamFileIn);
    readRecord(read2, bamFileIn);

    int r1_begin = read1.beginPos; // coordinates are 0-based!
    int r2_begin = read2.beginPos; // coordinates are 0-based!
    int r1_len = getAlignmentLengthInRef(read1);
    int r2_len = getAlignmentLengthInRef(read2);
    string r1_ref = toCString(contigNames(bamContextIn)[read1.rID]);
    string r2_ref = toCString(contigNames(bamContextIn)[read2.rID]);

    // check custom tags for read pair
    BamTagsDict tags(read1.tags);
    short phase = 0;
    int copy = 0; // default copy is 0 (base chrom.)
    // extract or set phase information
    unsigned idx_tag_xp = 0;
    if (findTagKey(idx_tag_xp, tags, "XP")) { // get phase info
      extractTagValue(phase, tags, idx_tag_xp);
    }
    else { // add custom field
      phase = rphase(); // pick random phase for read pair
      CharString tagXP = str(boost::format("XP:A:%s") % phase);
      appendTagsSamToBam(read1.tags, tagXP);
      appendTagsSamToBam(read2.tags, tagXP);
    }
    // extract or set copy information
    unsigned idx_tag_xc = 0;
    if (findTagKey(idx_tag_xc, tags, "XC")) { // get copy info
      extractTagValue(copy, tags, idx_tag_xc);
    }
    else { // add custom field
      CharString tagXC = str(boost::format("XC:A:%s") % copy);
      appendTagsSamToBam(read1.tags, tagXC);
      appendTagsSamToBam(read2.tags, tagXC);
    }

    // mutate reads
    // variants are accessible by chromosome name which should speed up lookup

    // find first candidate variant
    auto it_chr = variants.map_chr2pos2var.find(r1_ref);
    if (it_chr == variants.map_chr2pos2var.end()) {
      continue; // no variants for this reference sequence
    }
    auto it_pos = it_chr->second.lower_bound(r1_begin);

    // identify variants affecting read pair
    bool is_mutated = false;
    while (it_pos!=it_chr->second.end() && it_pos->first<r2_begin+r2_len) {
      Variant var = it_pos->second;
      is_mutated = (var.chr_copy == phase); // does read pair have same phase as variant?

      int r1_var_pos = var.pos-1 - r1_begin;
      if (r1_var_pos >= 0 && r1_var_pos < r1_len) { // read1 overlaps with variant
        map_var_cvg[var.id]++;
        if (is_mutated) { // mutate read1
          //cerr << str(boost::format("%s:%d\t%s->%s\n") % toCString(read1.qName) % r1_var_pos % var.alleles[0].c_str() % var.alleles[1].c_str());
          read1.seq[r1_var_pos] = var.alleles[1][0];
          map_var_vaf[var.id]++;
        }
      }
      int r2_var_pos = var.pos - r2_begin;
      if (r2_var_pos >= 0 && r2_var_pos < r2_len) { // read2 overlaps with variant
        map_var_cvg[var.id]++;
        if (is_mutated) { // mutate read2
          //cerr << str(boost::format("%s:%d\t%s->%s\n") % toCString(read2.qName) % r2_var_pos % var.alleles[0].c_str() % var.alleles[1].c_str());
          read2.seq[r2_var_pos] = var.alleles[1][0];
          map_var_vaf[var.id]++;
        }
      }
      ++it_pos;
    }

    // BAM output
    writeRecord(bamFileOut, read1);
    writeRecord(bamFileOut, read2);
  }

  fs_bamout.close();

  // write coverage and variant read counts for all variants to file
  string fn_varout = str(boost::format("%s.vars.csv") % fn_sam_out.c_str());
  ofstream fs_varout;
  fs_varout.open(fn_varout.c_str());
  for (Variant v : variants.vec_variants) {
    string id = str(boost::format("%s_%d:%ld") % v.chr % v.chr_copy % v.pos);
    unsigned alt = map_var_vaf[v.id];
    unsigned ref = map_var_cvg[v.id] - alt;
    string line = str(boost::format("%d\t%s\t%s\t%s\n") % v.idx_mutation % id % ref % alt);
    fs_varout << line;
  }
  fs_varout.close();
}

/** Spike in subclonal mutations to SAM input. */
// TODO: Ploidy could be clone-dependent.
void mutateReads(
  string fn_fq_out,
  string fn_sam_out,
  string fn_sam_in,
  vario::VariantSet &variants,
  vector<shared_ptr<Clone>> vec_vis_clones,
  map<string, vector<bool>> mm,
  vector<double> weights,
  string id_sample,
  const short ploidy,
  RandomNumberGenerator<> &rng,
  bool do_write_fastq)
{
  //sort clones by label
  size_t num_clones = vec_vis_clones.size();
  sort(vec_vis_clones.begin(), vec_vis_clones.end(),
       [](shared_ptr<Clone> a, shared_ptr<Clone> b) {
         return a->label < b->label;
       });
  vector<int> vec_clone_idx;
  vector<string> vec_clone_lbl;
  vector<double> vec_clone_weight;

  for (auto c : vec_vis_clones) {
    vec_clone_idx.push_back(c->index);
    vec_clone_lbl.push_back(c->label);
    //vec_clone_weight.push_back(c->weight); // is a parameter now -> multiple samples
  }
  // calculate cellular prevalence values (root: 1-sum)
  double s = 0.0;
  for (double w : weights) {
    vec_clone_weight.push_back(w);
    s += w;
  }
  //vec_clone_weight.insert(vec_clone_weight.begin(), 1.0-s);
  vec_clone_weight.push_back(1.0-s);

  // TODO: check verbosity setting
  fprintf(stdout, "-- %s --\n", id_sample.c_str());
  fprintf(stdout, "Baseline reads: %s\n", fn_sam_in.c_str());
  fprintf(stdout, "Creating bulk sample with the following prevalences:\n");
  for (int i=0; i<vec_clone_lbl.size(); i++) {
    fprintf(stdout, "clone %d: \"%s\"\t(%.4f)\n", vec_clone_idx[i], vec_clone_lbl[i].c_str(), vec_clone_weight[i]);
  }

  // store coverage depth and variant allele frequency for each variant
  map<string, unsigned> map_var_cvg;
  map<string, unsigned> map_var_vaf;

  // initialize random function
  function<short()> rphase = rng.getRandomFunctionInt(short(0), short(ploidy-1));
  // store coverage depth and variant allele frequency for each variant
  //vector<Variant> var_sorted = Variant::sortByPositionPoly(variants);
  for (Variant &var : variants.vec_variants) {
    map_var_cvg[var.id] = 0;
    map_var_vaf[var.id] = 0;
  }

  /*------------------*
   * process BAM file *
   *------------------*/

  function<int()> rand_idx = rng.getRandomIndexWeighted(vec_clone_weight);
  CharString bamFileName = getAbsolutePath(fn_sam_in.c_str());
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
  fs_bamout.open(fn_sam_out.c_str());
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
      /* healthy genome is generated from haploid reference now...remove this?
      if (pos_delim == string::npos) {
        fprintf(stderr, "[WARN] invalid sequence name: %s (missing copy info).\n", seq_name.c_str());
        appendValue(refNames, seqName);
        appendValue(refLengths, atoi(toCString(seqLen)));
        appendValue(headerOut, headRec);
      }
      else {*/
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
      //}
    }
    else {
      appendValue(headerOut, headRec);
    }
  }
  // add read group for each clone
  addCloneReadGroups(headerOut, id_sample, vec_clone_lbl);
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
  if (do_write_fastq) {
    fs_fq.open(fn_fq_out);
  }
  fs_log.open("bamio_bulk.log");
  BamAlignmentRecord read1, read2;
  while (!atEnd(bamFileIn)) {
    readRecord(read1, bamFileIn);
    readRecord(read2, bamFileIn);

    // pick random clone (weighted) to which next read pair belongs
    int r_idx = rand_idx();
    int c_idx = vec_clone_idx[r_idx];
    string c_lbl = vec_clone_lbl[c_idx];
    string r1_name = toCString(read1.qName);
    int r1_begin = read1.beginPos;
    int r2_begin = read2.beginPos;
    int r1_len = getAlignmentLengthInRef(read1);
    int r2_len = getAlignmentLengthInRef(read2);
    int r1_end = r1_begin + r1_len;
    int r2_end = r2_begin + r2_len;
    string r1_ref = toCString(contigNames(bamContextIn)[read1.rID]);
    string r2_ref = toCString(contigNames(bamContextIn)[read2.rID]);
    char r1_rc = seqan::hasFlagRC(read1) ? '+' : '-';
    char r2_rc = seqan::hasFlagRC(read2) ? '+' : '-';
    string r1_qual = toCString(read1.qual);
    string r2_qual = toCString(read2.qual);
    auto min_pos_read = (r1_begin < r2_begin ? r1_begin : r2_begin);
    auto max_pos_read = (r1_end > r2_end ? r1_end : r2_end);

    // check custom tags for read pair
    BamTagsDict tags(read1.tags);
    short phase = 0;
    int copy = 0; // default copy is 0 (base chrom.)
    // extract or set phase information
    unsigned idx_tag_xp = 0;
    if (findTagKey(idx_tag_xp, tags, "XP")) { // get phase info
      char char_phase;
      extractTagValue(char_phase, tags, idx_tag_xp);
      phase = char_phase - '0';
    }
    else { // add custom field
      phase = rphase(); // pick random phase for read pair
      CharString tagXP = str(boost::format("XP:A:%s") % phase);
      appendTagsSamToBam(read1.tags, tagXP);
      appendTagsSamToBam(read2.tags, tagXP);
    }
    // extract or set copy information
    unsigned idx_tag_xc = 0;
    if (findTagKey(idx_tag_xc, tags, "XC")) { // get copy info
      char char_copy;
      extractTagValue(char_copy, tags, idx_tag_xc);
      copy = char_copy - '0';
    }
    else { // add custom field
      string tagXC = str(boost::format("XC:A:%s") % copy);
      appendTagsSamToBam(read1.tags, tagXC);
      appendTagsSamToBam(read2.tags, tagXC);
    }
    fs_log << str(boost::format("%s:%s(%d,%d) -- %s:%s(%d,%d)\n")
      % r1_ref % r1_rc % r1_begin % (r1_begin+r1_len)
      % r2_ref % r2_rc % r2_begin % (r2_begin+r2_len));

    // modify read pair to match assigned clone, phase and copy
    //read1.qName += str(boost::format("-%d") % c_idx);
    //read2.qName += str(boost::format("-%d") % c_idx);
    CharString tagRG = str(boost::format("RG:Z:%s") % vec_clone_lbl[r_idx]);
    appendTagsSamToBam(read1.tags, tagRG);
    appendTagsSamToBam(read2.tags, tagRG);

    // mutate reads
    // variants are accessible by chromosome name which should speed up lookup

    // find first candidate variant
    auto it_chr = variants.map_chr2pos2var.find(r1_ref);
    if (it_chr == variants.map_chr2pos2var.end()) {
      continue; // no variants for this reference sequence
    }
    // TODO: misses mutations when occurring at same pos
    auto it_pos = it_chr->second.lower_bound(min_pos_read);

    // identify variants affecting read pair
    bool is_mutated = false;
    while (it_pos!=it_chr->second.end() && it_pos->first<=max_pos_read) {
      Variant var = it_pos->second;
      is_mutated = mm[c_lbl][var.idx_mutation]; // does read pair's clone carry mutation?
      is_mutated = is_mutated && (var.chr_copy == phase); // does variant phase match read pair phase?

      int r1_var_pos = var.pos - r1_begin;
      if (r1_var_pos >= 0 && r1_var_pos < r1_len) { // read1 overlaps with variant
        map_var_cvg[var.id]++;
        if (is_mutated) { // mutate read1
          map_var_vaf[var.id]++;
          fs_log << str(boost::format("%s:%d\t%s->%s\n") % toCString(read1.qName) % r1_var_pos % var.alleles[0].c_str() % var.alleles[1].c_str());
          read1.seq[r1_var_pos] = var.alleles[1][0];
        }
      }
      int r2_var_pos = var.pos+1 - r2_begin;
      if (r2_var_pos >= 0 && r2_var_pos < r2_len) { // read2 overlaps with variant
        map_var_cvg[var.id]++;
        if (is_mutated) { // mutate read2
          map_var_vaf[var.id]++;
          fs_log << str(boost::format("%s:%d\t%s->%s\n") % toCString(read2.qName) % r2_var_pos % var.alleles[0].c_str() % var.alleles[1].c_str());
          read2.seq[r2_var_pos] = var.alleles[1][0];
        }
      }

      ++it_pos;
    }

    // BAM output

    //CharString r1_ref_hap = refNames[read1.rID];
    //CharString r2_ref_hap = refNames[read2.rID];
    write(bamFileOut.iter, read1, bamContextMap, bamFileOut.format);
    write(bamFileOut.iter, read2, bamContextMap, bamFileOut.format);

    // FASTQ output
    if (do_write_fastq) {
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
      fs_fq << str(boost::format("@%s\n%s\n+\n%s\n") % read1.qName % read1.seq % r1_qual);
      fs_fq << str(boost::format("@%s\n%s\n+\n%s\n") % read2.qName % read2.seq % r2_qual);
    }
  }

  fs_bamout.close();
  if (do_write_fastq) {
    fs_fq.close();
  }
  fs_log.close();

  // write coverage and variant read counts for all variants to file
  string fn_varout = str(boost::format("%s.vars.csv") % fn_sam_out.c_str());
  ofstream fs_varout;
  fs_varout.open(fn_varout.c_str());
  for (Variant v : variants.vec_variants) {
    string id = str(boost::format("%s_%d:%ld") % v.chr % v.chr_copy % v.pos);
    unsigned alt = map_var_vaf[v.id];
    unsigned ref = map_var_cvg[v.id] - alt;
    string line = str(boost::format("%d\t%s\t%s\t%s\n") % v.id % id % ref % alt);
    fs_varout << line;
  }
  fs_varout.close();
}

/** DEPRECATED Spike in subclonal mutations to SAM input. */
// TODO: Ploidy could be clone-dependent.
void mutateReads(
  string fn_fq_out,
  string fn_sam_out,
  string fn_sam_in,
  vario::VariantSet &variants,
  treeio::Tree<Clone> &tree,
  vector<double> weights,
  string id_sample,
  const short ploidy,
  RandomNumberGenerator<> &rng,
  bool do_write_fastq)
{
  // extract info about visible clones
  vector<shared_ptr<Clone>> vec_vis_clones = tree.getVisibleNodes();
  size_t num_clones = vec_vis_clones.size();
  sort(vec_vis_clones.begin(), vec_vis_clones.end(),
       [](shared_ptr<Clone> a, shared_ptr<Clone> b) {
         return a->label < b->label;
       });
  vector<int> vec_clone_idx;
  vector<string> vec_clone_lbl;
  vector<double> vec_clone_weight;

  for (auto c : vec_vis_clones) {
    vec_clone_idx.push_back(c->index);
    vec_clone_lbl.push_back(c->label);
    //vec_clone_weight.push_back(c->weight); // is a parameter now -> multiple samples
  }
  // calculate cellular prevalence values (root: 1-sum)
  double s = 0.0;
  for (double w : weights) {
    vec_clone_weight.push_back(w);
    s += w;
  }
  //vec_clone_weight.insert(vec_clone_weight.begin(), 1.0-s);
  vec_clone_weight.push_back(1.0-s);

  // TODO: check verbosity setting
  fprintf(stdout, "-- %s --\n", id_sample.c_str());
  fprintf(stdout, "Creating bulk sample with the following prevalences:\n");
  for (int i=0; i<vec_clone_lbl.size(); i++) {
    fprintf(stdout, "clone %d: \"%s\"\t(%.4f)\n", vec_clone_idx[i], vec_clone_lbl[i].c_str(), vec_clone_weight[i]);
  }

  // get mutation matrix
  int num_nodes = tree.m_numNodes;
  vector<vector<bool>> mm(num_nodes, vector<bool>(variants.num_variants, false));
  tree.m_root->populateMutationMatrixRec(mm);

  // store coverage depth and variant allele frequency for each variant
  map<string, unsigned> map_var_cvg;
  map<string, unsigned> map_var_vaf;

  // initialize random function
  function<short()> rphase = rng.getRandomFunctionInt(short(0), short(ploidy-1));
  // store coverage depth and variant allele frequency for each variant
  //vector<Variant> var_sorted = Variant::sortByPositionPoly(variants);
  for (Variant &var : variants.vec_variants) {
    map_var_cvg[var.id] = 0;
    map_var_vaf[var.id] = 0;
  }

  /*------------------*
   * process BAM file *
   *------------------*/

  function<int()> rand_idx = rng.getRandomIndexWeighted(vec_clone_weight);
  CharString bamFileName = getAbsolutePath(fn_sam_in.c_str());
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
  fs_bamout.open(fn_sam_out.c_str());
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
      /* healthy genome is generated from haploid reference now...remove this?
      if (pos_delim == string::npos) {
        fprintf(stderr, "[WARN] invalid sequence name: %s (missing copy info).\n", seq_name.c_str());
        appendValue(refNames, seqName);
        appendValue(refLengths, atoi(toCString(seqLen)));
        appendValue(headerOut, headRec);
      }
      else {*/
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
      //}
    }
    else {
      appendValue(headerOut, headRec);
    }
  }
  // add read group for each clone
  addCloneReadGroups(headerOut, id_sample, vec_clone_lbl);
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
  if (do_write_fastq) {
    fs_fq.open(fn_fq_out);
  }
  fs_log.open("bamio_bulk.log");
  BamAlignmentRecord read1, read2;
  while (!atEnd(bamFileIn)) {
    readRecord(read1, bamFileIn);
    readRecord(read2, bamFileIn);

    // pick random clone (weighted) to which next read pair belongs
    int r_idx = rand_idx();
    int c_idx = vec_clone_idx[r_idx];
    int r1_begin = read1.beginPos;
    int r2_begin = read2.beginPos;
    int r1_len = getAlignmentLengthInRef(read1);
    int r2_len = getAlignmentLengthInRef(read2);
    string r1_ref = toCString(contigNames(bamContextIn)[read1.rID]);
    string r2_ref = toCString(contigNames(bamContextIn)[read2.rID]);
    char r1_rc = seqan::hasFlagRC(read1) ? '+' : '-';
    char r2_rc = seqan::hasFlagRC(read2) ? '+' : '-';
    string r1_qual = toCString(read1.qual);
    string r2_qual = toCString(read2.qual);

    // check custom tags for read pair
    BamTagsDict tags(read1.tags);
    short phase = 0;
    int copy = 0; // default copy is 0 (base chrom.)
    // extract or set phase information
    unsigned idx_tag_xp = 0;
    if (findTagKey(idx_tag_xp, tags, "XP")) { // get phase info
      char char_phase;
      extractTagValue(char_phase, tags, idx_tag_xp);
      phase = char_phase - '0';
    }
    else { // add custom field
      phase = rphase(); // pick random phase for read pair
      CharString tagXP = str(boost::format("XP:A:%s") % phase);
      appendTagsSamToBam(read1.tags, tagXP);
      appendTagsSamToBam(read2.tags, tagXP);
    }
    // extract or set copy information
    unsigned idx_tag_xc = 0;
    if (findTagKey(idx_tag_xc, tags, "XC")) { // get copy info
      char char_copy;
      extractTagValue(char_copy, tags, idx_tag_xc);
      copy = char_copy - '0';
    }
    else { // add custom field
      string tagXC = str(boost::format("XC:A:%s") % copy);
      appendTagsSamToBam(read1.tags, tagXC);
      appendTagsSamToBam(read2.tags, tagXC);
    }
    fs_log << str(boost::format("%s:%s(%d,%d) -- %s:%s(%d,%d)\n")
      % r1_ref % r1_rc % r1_begin % (r1_begin+r1_len)
      % r2_ref % r2_rc % r2_begin % (r2_begin+r2_len));

    // modify read pair to match assigned clone, phase and copy
    //read1.qName += str(boost::format("-%d") % c_idx);
    //read2.qName += str(boost::format("-%d") % c_idx);
    CharString tagRG = str(boost::format("RG:Z:%s") % vec_clone_lbl[r_idx]);
    appendTagsSamToBam(read1.tags, tagRG);
    appendTagsSamToBam(read2.tags, tagRG);

    // mutate reads
    // variants are accessible by chromosome name which should speed up lookup

    // find first candidate variant
    auto it_chr = variants.map_chr2pos2var.find(r1_ref);
    if (it_chr == variants.map_chr2pos2var.end()) {
      continue; // no variants for this reference sequence
    }
    auto it_pos = it_chr->second.lower_bound(r1_begin);

    // identify variants affecting read pair
    bool is_mutated = false;
    while (it_pos!=it_chr->second.end() && it_pos->first<r2_begin+r2_len) {
      Variant var = it_pos->second;
      is_mutated = mm[c_idx][var.idx_mutation]; // does read pair's clone carry mutation?
      is_mutated = is_mutated && (var.chr_copy == phase); // does variant phase match read pair phase?

      int r1_var_pos = var.pos-1 - r1_begin;
      if (r1_var_pos >= 0 && r1_var_pos < r1_len) { // read1 overlaps with variant
if (var.id == "m27") {
  fprintf(stderr, "%s\n", toCString(read1.qName));
}
        map_var_cvg[var.id]++;
        if (is_mutated) { // mutate read1
          map_var_vaf[var.id]++;
          fs_log << str(boost::format("%s:%d\t%s->%s\n") % toCString(read1.qName) % r1_var_pos % var.alleles[0].c_str() % var.alleles[1].c_str());
          read1.seq[r1_var_pos] = var.alleles[1][0];
        }
      }
      int r2_var_pos = var.pos - r2_begin;
      if (r2_var_pos >= 0 && r2_var_pos < r2_len) { // read2 overlaps with variant
if (var.id == "m27") {
  fprintf(stderr, "%s\n", toCString(read2.qName));
}
        map_var_cvg[var.id]++;
        if (is_mutated) { // mutate read2
          map_var_vaf[var.id]++;
          fs_log << str(boost::format("%s:%d\t%s->%s\n") % toCString(read2.qName) % r2_var_pos % var.alleles[0].c_str() % var.alleles[1].c_str());
          read2.seq[r2_var_pos] = var.alleles[1][0];
        }
      }

      ++it_pos;
    }

    // BAM output

    //CharString r1_ref_hap = refNames[read1.rID];
    //CharString r2_ref_hap = refNames[read2.rID];
    write(bamFileOut.iter, read1, bamContextMap, bamFileOut.format);
    write(bamFileOut.iter, read2, bamContextMap, bamFileOut.format);

    // FASTQ output
    if (do_write_fastq) {
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
      fs_fq << str(boost::format("@%s\n%s\n+\n%s\n") % read1.qName % read1.seq % r1_qual);
      fs_fq << str(boost::format("@%s\n%s\n+\n%s\n") % read2.qName % read2.seq % r2_qual);
    }
  }

  fs_bamout.close();
  if (do_write_fastq) {
    fs_fq.close();
  }
  fs_log.close();

  // write coverage and variant read counts for all variants to file
  string fn_varout = str(boost::format("%s.vars.csv") % fn_sam_out.c_str());
  ofstream fs_varout;
  fs_varout.open(fn_varout.c_str());
  for (Variant v : variants.vec_variants) {
    string id = str(boost::format("%s_%d:%ld") % v.chr % v.chr_copy % v.pos);
    unsigned alt = map_var_vaf[v.id];
    unsigned ref = map_var_cvg[v.id] - alt;
    string line = str(boost::format("%d\t%s\t%s\t%s\n") % v.id % id % ref % alt);
    fs_varout << line;
  }
  fs_varout.close();
}

void addCloneReadGroups(BamHeader &header, string id_sample, const vector<string> &vec_lbl) {
  BamHeaderRecord record;
  for (auto lbl : vec_lbl) {
    clear(record);
    record.type = seqan::BAM_HEADER_READ_GROUP;
    appendValue(record.tags, Pair<CharString>());
    assign(back(record.tags).i1, "ID", Exact());
    assign(back(record.tags).i2, lbl, Exact());
    appendValue(record.tags, Pair<CharString>());
    assign(back(record.tags).i1, "SM", Exact());
    assign(back(record.tags).i2, str(boost::format("%s_%s") % id_sample % lbl), Exact());
    appendValue(record.tags, Pair<CharString>());
    assign(back(record.tags).i1, "LB", Exact());
    assign(back(record.tags).i2, id_sample, Exact());
    appendValue(record.tags, Pair<CharString>());
    assign(back(record.tags).i1, "PL", Exact());
    assign(back(record.tags).i2, "Illumina", Exact());
    appendValue(record.tags, Pair<CharString>());
    assign(back(record.tags).i1, "PU", Exact());
    assign(back(record.tags).i2, "HighSeq2500", Exact());

    appendValue(header, record);
  }
}

} // namespace bamio
