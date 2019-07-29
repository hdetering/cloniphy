#include "../seqio.hpp"
#include "ChromosomeInstance.hpp"

using namespace std;
using boost::icl::interval;

namespace seqio {

ChromosomeInstance::ChromosomeInstance() : length(0) {}

ChromosomeInstance::ChromosomeInstance (
  const ChromosomeReference ref,
  const char gl_allele
)
: length(ref.length) {
  // inititally chromosome consists of a single SegmentCopy
  SegmentCopy seg_copy(0, this->length, gl_allele);
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

vector<seg_mod_t> ChromosomeInstance::amplifyRegion (
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
  //bool is_right_bkp = false; // has the end of insertion been reached?
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
  // make sure insert length is exact
  len_bp = (bkp_end - bkp_start);

  // identify first affected SegmentCopy
  unsigned long pos_bp = 0; // current position
  unsigned long seg_len = 0; // length of current SegmentCopy
  char seg_allele = 'A'; // reference source allele of current SegmentCopy
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
    seg_allele = it_seg->gl_allele;
    //is_right_bkp = (bkp_end <= pos_bp+seg_len);
    // create new SegmentCopy
    unsigned long seg_new_start = it_seg->ref_start + (bkp_start - pos_bp);
    unsigned long seg_new_end = pos_bp+seg_len <= bkp_end ? it_seg->ref_end : seg_new_start+len_bp;
    assert( seg_new_start < seg_new_end );
    SegmentCopy seg_new(seg_new_start, seg_new_end, seg_allele);
    lst_seg_new.push_back(seg_new);
    vec_seg_mod.push_back(make_tuple(seg_new.id, it_seg->id, seg_new_start, seg_new_end));

    // will new SegmentCopies be inserted at first breakpoint?
    if (!is_forward ) {
      if ( bkp_start > pos_bp) { // split this SegmentCopy at breakpoint
        // head of SegmentCopy
        TCoord head_len = bkp_start - pos_bp;
        TCoord head_start = it_seg->ref_start;
        TCoord head_end = it_seg->ref_start + head_len;
        SegmentCopy seg_head(head_start, head_end, seg_allele);
        vec_seg_mod.push_back(make_tuple(seg_head.id, it_seg->id, head_start, head_end));
        // tail of SegmentCopy
        TCoord tail_start = head_end;
        TCoord tail_end = it_seg->ref_end;
        SegmentCopy seg_tail(tail_start, tail_end, seg_allele);
        vec_seg_mod.push_back(make_tuple(seg_tail.id, it_seg->id, tail_start, tail_end));

        // replace existing SegmentCopy
        it_seg = this->lst_segments.erase(it_seg);
        this->lst_segments.insert(it_seg, seg_head);
        it_seg_bkp_left = this->lst_segments.insert(it_seg, seg_tail);
        // continue processing with new tail SegmentCopy
        it_seg = it_seg_bkp_left;
        pos_bp += head_len;
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
      SegmentCopy seg_copy(it_seg->ref_start, it_seg->ref_end, it_seg->gl_allele);
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
      TCoord seg_new_start = it_seg->ref_start;
      TCoord seg_new_end = seg_new_start + (bkp_end-pos_bp);
      char seg_new_allele = it_seg->gl_allele;
      assert( seg_new_start < seg_new_end );
      SegmentCopy seg_new(seg_new_start, seg_new_end, seg_new_allele);
      lst_seg_new.push_back(seg_new);
      vec_seg_mod.push_back(make_tuple(seg_new.id, it_seg->id, seg_new_start, seg_new_end));
    }
    if (is_forward ) {
      if (bkp_end < pos_bp+seg_len) { // split this SegmentCopy at breakpoint
        // head of SegmentCopy
        TCoord head_len = bkp_end - pos_bp;
        TCoord head_start = it_seg->ref_start;
        TCoord head_end = it_seg->ref_start + head_len;
        SegmentCopy seg_head(head_start, head_end, it_seg->gl_allele);
        vec_seg_mod.push_back(make_tuple(seg_head.id, it_seg->id, head_start, head_end));

        // tail of SegmentCopy
        TCoord tail_start = head_end;
        TCoord tail_end = it_seg->ref_end;
        SegmentCopy seg_tail(tail_start, tail_end, it_seg->gl_allele);
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
  TCoord start_bp = start_rel * this->length;
  TCoord len_bp = len_rel * this->length;
  // physical coordinates
  TCoord bkp_start, bkp_end;
  if (is_forward) {
    bkp_start = start_bp;
    bkp_end = is_telomeric ? this->length : start_bp+len_bp;
  } else { // !is_forward
    bkp_end = start_bp;
    bkp_start = is_telomeric ? 0 : bkp_end-len_bp;
  }
  // make sure insert length is exact
  len_bp = (bkp_end - bkp_start);

  // identify first affected SegmentCopy
  TCoord pos_bp = 0; // current position
  TCoord seg_len = 0; // length of current SegmentCopy
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
    TCoord head_len = bkp_start - pos_bp;
    TCoord head_start = it_seg->ref_start;
    TCoord head_end = it_seg->ref_start + head_len;
    assert( head_start < head_end );
    SegmentCopy seg_head(head_start, head_end, it_seg->gl_allele);
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
    TCoord tail_len = pos_bp + seg_len - bkp_end;
    TCoord tail_start = it_seg->ref_end - tail_len;
    TCoord tail_end = it_seg->ref_end;
    assert( tail_start < tail_end );
    SegmentCopy seg_tail(tail_start, tail_end, it_seg->gl_allele);
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
    SegmentCopy seg_new(seg_old.ref_start, seg_old.ref_end, seg_old.gl_allele);
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

bool
ChromosomeInstance::indexSegmentCopies (
  TSegMap& imap_segments
) const
{
  typedef set<SegmentCopy> TSegSet;
  for (const SegmentCopy seg : lst_segments) {
    TCoord p1 = seg.ref_start;
    TCoord p2 = seg.ref_end;
    imap_segments += make_pair(interval<TCoord>::right_open(p1,p2), TSegSet({seg}));
  }
}

} // namespace seqio