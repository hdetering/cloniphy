#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include <boost/format.hpp>
using boost::format;
using boost::str;
#include <fstream>
#include "../core/clone.cpp"
#include "../core/random.hpp"
#include "../core/treeio.hpp"
using namespace treeio;

struct FixtureTreeio {
  FixtureTreeio() {
    BOOST_TEST_MESSAGE( "setup fixure" );
  }
  ~FixtureTreeio() {
    BOOST_TEST_MESSAGE( "teardown fixture" );
  }

  long seed = 123456789;
  RandomNumberGenerator rng = RandomNumberGenerator(seed);
};

BOOST_FIXTURE_TEST_SUITE( treeio, FixtureTreeio )

/* generate random tree topology */
BOOST_AUTO_TEST_CASE( random_clone_tree )
{
  boost::timer::auto_cpu_timer t;
  function<double()> random_dbl = rng.getRandomFunctionReal(0.0, 1.0);
  function<double()> random_gamma = rng.getRandomFunctionGamma(2.0, 0.25);
  int num_clones = 5;
  Tree<Clone> tree(num_clones);
  tree.generateRandomTopologyInternalNodes(random_dbl);
  tree.varyBranchLengths(random_gamma);
  int num_mutations = 1010;
  int num_transmuts =   10;
  tree.dropSomaticMutations(num_mutations, num_transmuts, rng);

  BOOST_CHECK( tree.m_numVisibleNodes == num_clones+1 );
  BOOST_CHECK( tree.m_numNodes == tree.m_vecNodes.size() );

  // assign random weights
  vector<double> w = rng.getRandomProbs(num_clones);
  double c = 0.1; // contamination with normal cells
  for (auto &p : w) p*=(1-c);
  w.insert(w.begin(), c);
  double sum_probs = 0.0;
  for (auto p : w) sum_probs += p;
  BOOST_CHECK( sum_probs == 1.0 );
  BOOST_TEST_MESSAGE( str(format("sum_probs: %0.4f") % sum_probs) );
  tree.assignWeights(w);

  BOOST_TEST_MESSAGE( "Writing resulting tree to files:" );
  BOOST_TEST_MESSAGE( "  random_clone_tree.dot" );
  ofstream fs_dot;
  fs_dot.open("random_clone_tree.dot");
  tree.printDot(tree.m_root, fs_dot);
  fs_dot.close();
  BOOST_TEST_MESSAGE( "  random_clone_tree.tre" );
  ofstream fs_nwk;
  fs_nwk.open("random_clone_tree.tre");
  tree.printNewick(fs_nwk);
  fs_nwk.close();
}

/* generate random tree topology */
BOOST_AUTO_TEST_CASE( random_sample_tree )
{
  boost::timer::auto_cpu_timer t;
  function<double()> random_dbl = rng.getRandomFunctionReal(0.0, 1.0);
  function<double()> random_gamma = rng.getRandomFunctionGamma(2.0, 0.5);
  int num_samples = 5;
  Tree<Clone> tree(num_samples);
  tree.generateRandomTopologyLeafsOnly(random_dbl);
  tree.varyBranchLengths(random_gamma);
  int num_mutations = 1010;
  int num_transmuts =   10;
  tree.dropSomaticMutations(num_mutations, num_transmuts, rng);

  BOOST_CHECK( tree.m_numVisibleNodes == num_samples+1 );
  BOOST_CHECK( tree.m_numNodes == tree.m_vecNodes.size() );

  // assign random weights
  double c = 0.1; // contamination with normal cells
  vector<double> w = { c, 0.4, 0.3, 0.1, 0.075,	0.025 };
  /*vector<double> w = rng.getRandomProbs(num_samples);
  w.insert(w.begin(), c);
  for_each(w.begin()+1, w.end(), [&](double &p) { p*=(1-c); });*/
  double sum_probs = 0.0;
  for (auto p : w) sum_probs += p;
  BOOST_CHECK( sum_probs == 1.0 );
  tree.assignWeights(w);

  BOOST_TEST_MESSAGE( "Writing resulting tree to output files:" );
  BOOST_TEST_MESSAGE( "  random_sample_tree.dot" );
  ofstream fs_dot;
  fs_dot.open("random_sample_tree.dot");
  tree.printDot(tree.m_root, fs_dot);
  fs_dot.close();
  BOOST_TEST_MESSAGE( "  random_sample_tree.tre" );
  ofstream fs_nwk;
  fs_nwk.open("random_sample_tree.tre");
  tree.printNewick(fs_nwk);
  fs_nwk.close();

  BOOST_TEST_MESSAGE( "Building mutation matrix from tree..." );
  int num_nodes = tree.m_numNodes;
  vector<vector<bool>> mm(num_nodes, vector<bool>(num_mutations, false));
  tree.m_root->populateMutationMatrixRec(mm);
}

/* read user-specified tree from file */
BOOST_AUTO_TEST_CASE ( read )
{
  string fn = "data/bulk/bulk.tre";

  BOOST_TEST_MESSAGE( "Reading clone tree from file:\n  " << fn << "\n");

  Tree<Clone> tree(fn);
  tree._printNodes();
  vector<shared_ptr<Clone>> clones = tree.getVisibleNodes();

  BOOST_TEST_MESSAGE( "\nVisible nodes (clones):" );

  int i = 0;
  for (auto c : clones) {
    BOOST_TEST_MESSAGE( str(format("  %d: '%s'") % i % c->label ) );
    i++;
  }
}

BOOST_AUTO_TEST_SUITE_END()
