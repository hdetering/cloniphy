#include "Variant.hpp"

using namespace std;

namespace vario {

Variant::Variant(std::string id, std::string chr, unsigned long pos)
 : id(id), chr(chr), pos(pos) {}

/* Variant *
 *---------*/

Variant::Variant ()
: id(""),
  chr(""),
  pos(0),
  alleles(0),
  idx_mutation(0),
  rel_pos(0.0),
  is_somatic(false),
  is_het(true),
  is_error(false)
{}
Variant::~Variant () {}

bool Variant::operator< (const Variant &other) const {
  return rel_pos < other.rel_pos;
}

vector<Variant> Variant::sortByPosition(const vector<Variant> &variants) {
  vector<Variant> variantsCopy = variants;
  sort(variantsCopy.begin(), variantsCopy.end());
  return variantsCopy;
}

vector<Variant> Variant::sortByPositionLex(const vector<Variant> &variants) {
  vector<Variant> variantsCopy = variants;
  sort(variantsCopy.begin(), variantsCopy.end(),
      [](const Variant &a, const Variant &b) -> bool {
        if (a.chr < b.chr) return true;
        if (a.chr == b.chr) return a.pos < b.pos;
        return false;
      });
  return variantsCopy;
}

vector<Variant> Variant::sortByPositionRef(const vector<Variant> &variants) {
  vector<Variant> variantsCopy = variants;
  sort(variantsCopy.begin(), variantsCopy.end(),
      [](const Variant &a, const Variant &b) -> bool {
        return a.pos < b.pos;
      });
  return variantsCopy;
}

bool Variant::isSnv() {
  for (vector<string>::iterator allele=this->alleles.begin(); allele!=this->alleles.end(); ++allele) {
    if ((*allele).size()>1) { return false; }
  }
  return true;
}


} // namespace vario