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
  return static_cast<float>(m_vec_mutations.size());
}

void Clone::replace(Clone *cloneToReplace) {
  // inherit properties
  this->m_vec_mutations = cloneToReplace->m_vec_mutations;
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

// returns filenames of genome files in clone2fn
void Clone::mutateGenome(const Genome &ref_genome,
                         const vector<Mutation> &mutations,
                         vector<Variant> &variants,
                         vector<vector<short> > &mutMatrix,
                         map<Clone*, string>& clone2fn) {
cerr << "\nGenerating genome for " << *this << ", " << this->m_vec_mutations.size() << " mutations" << endl;

  // collect mutations hierarcically ascending in the tree
  vector<Mutation> my_mutations;
  for (Clone *c=this; c->parent!=0; c=c->parent) {
    // inherit ancestral mutations
    for (vector<int>::iterator i=c->m_vec_mutations.begin(); i!=c->m_vec_mutations.end(); ++i) {
      my_mutations.push_back(mutations[*i]);
      mutMatrix[this->index][*i] = mutations[*i].copy+1; // remember genotype (easier VCF export)
    }
  }

  // assign clone id (either by label or index)
  string clone_id = "";
  if (this->label.size()>0)
    clone_id = this->label;
  else
    clone_id = boost::str(boost::format("clone_%02d") % this->index);

  // rename sequences adding clone id
  // (maybe not such a good idea afterall...)
  //for (unsigned j=0; j<clone_genome.records.size(); ++j) {
  //  clone_genome.records[j].id += '-' + clone_id;
  //}

  // write clone genome to file
  string filename = boost::str(boost::format("%s_genome.fa") % clone_id);
  ofstream outfile;
  outfile.open(filename.c_str());
  //seqio::writeFasta(my_genome.records, outfile);
  vario::applyVariantsStream(ref_genome, my_mutations, variants, outfile);
  outfile.close();
  clone2fn[this] = filename;
//cerr << endl;
}

void Clone::applyMutations(const vector<Mutation> &my_mutations,
                          SubstitutionModel model,
                          Genome &my_genome,
                          boost::function<double()>& rng)
{
  // sort mutations by position
  std::vector<Mutation> mut_sorted = Mutation::sortByPosition(my_mutations);
  // apply mutations in order of reference position
  for (std::vector<Mutation>::iterator i=mut_sorted.begin(); i!=mut_sorted.end(); ++i) {
    Mutation m = *i;
    Variant var = *(new Variant());
    Genotype gt = *(new Genotype());
    m.apply(my_genome, model, rng, var, gt);
    //m_vec_variants.push_back(var); // append to own variants
    m_vec_genotypes.push_back(gt);
  }
}