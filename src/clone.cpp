#include "clone.hpp"
#include <ostream>

Clone::Clone() {
  this->label = -1;
  this->freq = 0.0;
  this->is_healthy = false;
  this->is_visible = false;
}

Clone::~Clone() {}

bool Clone::isLeaf() {
  return m_vecChildren.size() == 0;
}

std::vector<Clone *> Clone::getChildren() {
  return m_vecChildren;
}

/** Represents a clone node in a clonal tree */
float Clone::distanceToParent() {
  return static_cast<float>(m_vecMutations.size());
}

void Clone::replace(Clone *cloneToReplace) {
  // inherit properties
  this->m_vecMutations = cloneToReplace->m_vecMutations;
  // "adopt" children
  for (unsigned j=0; j<cloneToReplace->m_vecChildren.size(); j++) {
    if (cloneToReplace->m_vecChildren[j] != this) { // don't add child as its own child
      this->m_vecChildren.push_back(cloneToReplace->m_vecChildren[j]);
      cloneToReplace->m_vecChildren[j]->parent = this;
    }
  }
  // replace parent
  this->parent = cloneToReplace->parent;
  // replace clone in parent's children
  for (unsigned j=0; j<cloneToReplace->parent->m_vecChildren.size(); j++) {
    if (cloneToReplace->parent->m_vecChildren[j] == cloneToReplace) {
      cloneToReplace->parent->m_vecChildren[j] = this;
    }
  }
}

void Clone::mutateGenome(std::vector<SeqRecord> &genome, const std::vector<long>& cumStart, const std::vector<Mutation> &mutations) {
std::cerr << std::endl;
  // collect ancestral mutations
  std::vector<int> mut_ids;
  for (Clone *c=this; c->parent!=0; c=c->parent) {
    mut_ids.insert(mut_ids.end(), c->m_vecMutations.begin(), c->m_vecMutations.end());
  }
  std::vector<Mutation> myMutations;
  for (std::vector<int>::iterator i=mut_ids.begin(); i!=mut_ids.end(); ++i) {
    myMutations.push_back(mutations[*i]);
  }
  // sort mutations by position
  std::vector<Mutation> mutSorted = Mutation::sortByPosition(myMutations);
  // apply mutations in order of reference position
  for (std::vector<Mutation>::iterator i=mutSorted.begin(); i!=mutSorted.end(); ++i) {
    Mutation m = *i;
//std::cerr << (*i).absPos << ", " << (*i).offset << std::endl;
    unsigned s=0;
    while (s<cumStart.size() && s<m.absPos) { s++; }
    long loc_pos = m.absPos - cumStart[s-1];
    Nuc nuc_old = SeqIO::charToNuc(genome[s-1].seq[loc_pos]);
    Nuc nuc_new = static_cast<Nuc>((static_cast<int>(nuc_old) + m.offset) % 4); // TODO: this could be more generic (get rid of the hard-coded 4)
    genome[s-1].seq[loc_pos] = SeqIO::nucToChar(nuc_new);
std::cerr << "<Mutation(abs_pos=" << m.absPos << ",offset=" << m.offset << ")> mutating " << nuc_old << " to " << nuc_new << "" << std::endl;
  }
//std::cerr << "applying a total of " << mut_ids.size() << " mutations to <Clone(label=" << this->label << ")>" << std::endl;
}
