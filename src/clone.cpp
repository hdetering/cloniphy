#include "clone.hpp"
#include <ostream>

Clone::Clone() {
  this->label = -1;
  this->freq = 0.0;
  this->is_healthy = false;
  this->is_visible = false;
}

Clone::~Clone() {}

void Clone::setParent(Clone *parent) {
  this->parent = parent;
  parent->m_vecChildren.push_back(this);
}

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

void Clone::mutateGenome(std::vector<SeqRecord> &genome, const std::vector<unsigned long>& cumStart, const std::vector<Mutation> &mutations, std::vector<short> &mutMatrixRow) {
  // collect ancestral mutations
  std::vector<int> mut_ids;
  for (Clone *c=this; c->parent!=0; c=c->parent) {
    mut_ids.insert(mut_ids.end(), c->m_vecMutations.begin(), c->m_vecMutations.end());
  }
  std::vector<Mutation> myMutations;
  for (std::vector<int>::iterator i=mut_ids.begin(); i!=mut_ids.end(); ++i) {
    myMutations.push_back(mutations[*i]);
    mutMatrixRow[mutations[*i].id] = 1;
  }
  // sort mutations by position
  std::vector<Mutation> mutSorted = Mutation::sortByPosition(myMutations);
  // apply mutations in order of reference position
  unsigned numSeqs = cumStart.size();
  for (std::vector<Mutation>::iterator i=mutSorted.begin(); i!=mutSorted.end(); ++i) {
    Mutation m = *i;
//std::cerr << (*i).absPos << ", " << (*i).offset << std::endl;
    unsigned s=0;
    while (s<cumStart.size() && s<=m.absPos) { s++; }
    unsigned long loc_pos = m.absPos - cumStart[s-1];
    unsigned targetSeqIndex = (s-1)+(numSeqs*m.copy);
    char old_base = genome[targetSeqIndex].seq[loc_pos];
    char new_base = SeqIO::shiftNucleotide(old_base, m.offset);
    genome[targetSeqIndex].seq[loc_pos] = new_base;
std::cerr << "<Mutation(abs_pos=" << m.absPos << ",offset=" << m.offset << ",copy=" << m.copy << ")> mutating " << old_base << " to " << new_base << "" << std::endl;
  }
//std::cerr << "applying a total of " << mut_ids.size() << " mutations to <Clone(label=" << this->label << ")>" << std::endl;
std::cerr << std::endl;
}
