#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include <boost/function.hpp>
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
  RandomNumberGenerator<> rng = RandomNumberGenerator<>(seed);
};

BOOST_FIXTURE_TEST_SUITE( treeio, FixtureTreeio )

/* generate random tree topology */
BOOST_AUTO_TEST_CASE( random_clone_tree )
{
  boost::function<double()> random_dbl = rng.getRandomFunctionDouble(0.0, 1.0);
  int num_clones = 10;
  Tree<Clone> tree(num_clones);
  tree.generateRandomTopologyInternalNodes(random_dbl);

  BOOST_CHECK( tree.m_numNodes == num_clones );

  BOOST_TEST_MESSAGE( "Writing resulting tree to file 'random_clone_tree.dot'." );
  ofstream fs_dot;
  fs_dot.open("random_clone_tree.dot");
  tree.printDot(tree.m_root, fs_dot);
  fs_dot.close();
}

/* generate random tree topology */
BOOST_AUTO_TEST_CASE( random_sample_tree )
{
  boost::function<double()> random_dbl = rng.getRandomFunctionDouble(0.0, 1.0);
  int num_samples = 10;
  Tree<Clone> tree(num_samples);
  tree.generateRandomTopologyLeafsOnly(random_dbl);

  BOOST_CHECK( tree.m_numNodes == num_samples );

  BOOST_TEST_MESSAGE( "Writing resulting tree to file 'random_sample_tree.dot'." );
  ofstream fs_dot;
  fs_dot.open("random_sample_tree.dot");
  tree.printDot(tree.m_root, fs_dot);
  fs_dot.close();
}


BOOST_AUTO_TEST_SUITE_END()
