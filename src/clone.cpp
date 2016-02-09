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

// returns filenames of genome files in clone2fn
void Clone::mutateGenome(const Genome &ref_genome,
                         const vector<Mutation> &mutations,
                         SubstitutionModel model,
                         vector<Variant> &variants,
                         vector<vector<short> > &mutMatrix,
                         boost::function<float()>& rng,
                         map<Clone*, string>& clone2fn) {
cerr << "\nGenerating genome for " << *this << ", " << this->m_vecMutations.size() << " mutations" << endl;

  { // scope makes sure variables are deleted after use (reduces memory consumption)
    Genome my_genome = ref_genome; // TODO: check memory footprint
    // collect ancestral variants
    vector<Variant> anc_variants;
    vector<Genotype> anc_genotypes;
    if (this->parent != 0) {
      for (Clone *c=this->parent; c->parent!=0; c=c->parent) {
        // inherit ancestral mutations
        for (vector<int>::iterator i=c->m_vecMutations.begin(); i!=c->m_vecMutations.end(); ++i)
          mutMatrix[this->index][*i] = mutMatrix[c->index][*i];
        anc_variants.insert(anc_variants.end(), c->m_vec_variants.begin(), c->m_vec_variants.end());
        anc_genotypes.insert(anc_genotypes.end(), c->m_vec_genotypes.begin(), c->m_vec_genotypes.end());
      }
    }
    // apply ancestral variants
    vario::applyVariants(my_genome, anc_variants, anc_genotypes);

    // apply mutations for self, creating new variants
    vector<Mutation> my_mutations = vector<Mutation>(0);
    for (vector<int>::iterator i=this->m_vecMutations.begin(); i!=this->m_vecMutations.end(); ++i) {
      my_mutations.push_back(mutations[*i]);
      // remember the chromosome copy the mutation applies to
      mutMatrix[this->index][mutations[*i].id] = mutations[*i].copy+1;
    }
    this->applyMutations(my_mutations, model, my_genome, rng);

    // remember variants generated for this clone
    for (std::vector<Variant>::iterator var=m_vec_variants.begin(); var!=m_vec_variants.end(); ++var)
      variants.push_back(*var); // append to global variants

    // am I a clone whose genome is to be output?
    if (this->is_visible) {
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
      seqio::writeFasta(my_genome.records, outfile);
      outfile.close();
      clone2fn[this] = filename;
    }
  } // end of arbitrary scope

  // recurse for children
  for (unsigned i=0; i<this->m_vecChildren.size(); i++) {
    m_vecChildren[i]->mutateGenome(ref_genome, mutations, model, variants, mutMatrix, rng, clone2fn);
  }
//std::cerr << "applying a total of " << mut_ids.size() << " mutations to <Clone(label=" << this->label << ")>" << std::endl;
cerr << endl;
}

void Clone::applyMutations(const vector<Mutation> &my_mutations,
                          SubstitutionModel model,
                          Genome &my_genome,
                          boost::function<float()>& rng)
{
  // sort mutations by position
  std::vector<Mutation> mut_sorted = Mutation::sortByPosition(my_mutations);
  // apply mutations in order of reference position
  for (std::vector<Mutation>::iterator i=mut_sorted.begin(); i!=mut_sorted.end(); ++i) {
    Mutation m = *i;
    Variant var = *(new Variant());
    Genotype gt = *(new Genotype());
    m.apply(my_genome, model, rng, var, gt);
    m_vec_variants.push_back(var); // append to own variants
    m_vec_genotypes.push_back(gt);
  }
}
