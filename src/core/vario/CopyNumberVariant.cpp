#include "CopyNumberVariant.hpp"
#include <fstream>

using namespace std;

namespace vario {

CopyNumberVariant::CopyNumberVariant ()
: id(0),
  is_wgd(false),
  is_deletion(false),
  is_chr_wide(false),
  is_telomeric(false),
  is_first_arm(false),
  is_downstream(true),
  len_rel(0.0),
  start_rel(0.0),
  ref_pos_begin(0),
  ref_pos_end(0),
  len_abs(0),
  ref_chr("")
{}

ostream& operator<<(ostream& lhs, const CopyNumberVariant& cnv) {
  lhs << cnv.id;
  lhs << "\t" <<(cnv.is_wgd ? "WGD" : (cnv.is_deletion ? "DEL" : "CPY"));
  lhs << "\t" << (cnv.is_chr_wide ? "chr" : (cnv.is_telomeric ? "tel" : "foc"));
  lhs << "\t" << cnv.ref_chr;
  lhs << "\t" << cnv.ref_pos_begin + 1;
  lhs << "\t" << cnv.ref_pos_end + 1;
  lhs << "\t" << cnv.len_abs;
  lhs << "\t" << cnv.ref_allele;
  lhs << "\t" << (cnv.is_downstream ? "+" : "-");
  lhs << "\n";
  return lhs;
}

} // namespace vario