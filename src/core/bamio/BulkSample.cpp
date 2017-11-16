#include "BulkSample.hpp"

using namespace std;

namespace bamio {
  

BulkSample::BulkSample() {}

BulkSample::BulkSample(
  const string id_sample,
  const map<string, double> clone_weight
)
: id(id_sample),
  m_clone_weight(clone_weight)
{}

vector<string>
BulkSample::getCloneLabels ()
{
  vector<string> vec_clone_lbl;

  for (auto lbl_w : this->m_clone_weight) {
    vec_clone_lbl.push_back(lbl_w.first);
  }

  return vec_clone_lbl;
}


} // namespace bamio