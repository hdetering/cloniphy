#include "clone.hpp"
#include <boost/format.hpp>
#include <ostream>

using namespace std;

Clone::Clone() {
  this->label = -1;
  this->freq = 0.0;
  this->is_healthy = false;
  this->is_visible = false;
}

Clone::~Clone() {}

void Clone::setParent(shared_ptr<Clone> child, shared_ptr<Clone> parent) {
  child->parent = parent;
  parent->m_vecChildren.push_back(child);
}

bool Clone::isLeaf() {
  return m_vecChildren.size() == 0;
}

vector<shared_ptr<Clone>> Clone::getChildren() {
  return m_vecChildren;
}

void Clone::populateMutationMatrixRec(vector<vector<bool>> &mm) {
  // make sure the mutations matrix includes this clone
  if (mm.size() <= this->index ) {
    fprintf(stderr, "[ERROR] mutation matrix does not include index %d.\n", this->index);
    return;
  }
  // inherit mutations from parent
  if (this->parent)
    for (auto i=0; i<mm[this->index].size(); ++i)
      mm[this->index][i] = mm[this->parent->index][i];
  for (auto m : this->m_vec_mutations)
    mm[this->index][m] = true;
  for (auto child : this->m_vecChildren)
    child->populateMutationMatrixRec(mm);
}

void Clone::populateMutationOccRec(vector<vector<bool>> &mm) {
  // make sure the mutations matrix includes this clone
  if (mm.size() <= this->index ) {
    fprintf(stderr, "[ERROR] mutation matrix does not include index %d.\n", this->index);
    return;
  }
  for (auto m : this->m_vec_mutations)
    mm[this->index][m] = true;
  for (auto child : this->m_vecChildren)
    child->populateMutationMatrixRec(mm);
}

/** Represents a clone node in a clonal tree */
float Clone::distanceToParent() {
  return static_cast<float>(m_vec_mutations.size());
}

void Clone::replace(shared_ptr<Clone> cloneToReplace, shared_ptr<Clone> newClone) {
  // inherit properties
  newClone->m_vec_mutations = cloneToReplace->m_vec_mutations;
  // "adopt" children
  for (unsigned j=0; j<cloneToReplace->m_vecChildren.size(); j++) {
    if (cloneToReplace->m_vecChildren[j] != newClone) { // don't add child as its own child
      newClone->m_vecChildren.push_back(cloneToReplace->m_vecChildren[j]);
      cloneToReplace->m_vecChildren[j]->parent = newClone;
    }
  }
  // replace parent
  newClone->parent = cloneToReplace->parent;
  // replace clone in parent's children
  for (unsigned j=0; j<cloneToReplace->parent->m_vecChildren.size(); j++) {
    if (cloneToReplace->parent->m_vecChildren[j] == cloneToReplace) {
      cloneToReplace->parent->m_vecChildren[j] = newClone;
    }
  }
}
