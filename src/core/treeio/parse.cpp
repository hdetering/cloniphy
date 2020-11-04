#include "parse.hpp"
#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <fstream>
// #include <stdexcept>
// #include <streambuf>

using namespace std;

/*--------------------------------------*/
/*            parser logic              */
/*--------------------------------------*/

BOOST_FUSION_ADAPT_STRUCT(
  treeio::parse::node,
  (treeio::parse::children_vector, children)
  (string, label)
  (double, length)
)

namespace treeio {
namespace parse {

  namespace qi = boost::spirit::qi;
  namespace ascii = boost::spirit::ascii;

  /** Represents the Newick format grammar.
    *
    * For details see http://evolution.genetics.washington.edu/phylip/newick_doc.html
    * All nodes must have labels.
    */
  template <typename Iterator>
  struct newick_grammar : qi::grammar<Iterator, node(), ascii::space_type>
  {
    newick_grammar() : newick_grammar::base_type(tree)
    {
      // use %= to assign the result of the parse to the string
      label %= quoted_label | unquoted_label;
      // unquoted labels may not contain formatting characters
      unquoted_label = qi::lexeme[*(qi::char_ - ':' - ',' - ';' - '(' - ')' - '[' - ']' - '\'')];
      // quoted labels can contain any characters
      quoted_label = '\'' >> qi::lexeme[*(qi::char_ - '\'')] >> '\'';

      // use %= to assign the result of the parse to the double
      branch_length %= ':' >> qi::double_;

      // Parser return values will be assigned to corresponding struct elements
      // in the order specified!
      subtree = -descendant_list >> -label >> -branch_length;

      // Descendant list is a vector of node, we just push back the
      // created nodes into the vector
      descendant_list = '(' >> subtree >> *(',' >> subtree ) >> ')';

      tree %= subtree >> ';';

      // prepare parsers for debugging
      BOOST_SPIRIT_DEBUG_NODE(label);
      BOOST_SPIRIT_DEBUG_NODE(branch_length);
      BOOST_SPIRIT_DEBUG_NODE(subtree);
      BOOST_SPIRIT_DEBUG_NODE(descendant_list);
      BOOST_SPIRIT_DEBUG_NODE(tree);
      // activate debugging output for parsers
      //qi::debug(tree);
      //qi::debug(subtree);
      //qi::debug(descendant_list);
      //qi::debug(label);
      //qi::debug(branch_length);
    }

    private:
      /* grammar rules, typed by element they create */
      qi::rule<Iterator, node(), ascii::space_type> tree, subtree;
      qi::rule<Iterator, vector<node>(), ascii::space_type> descendant_list;
      qi::rule<Iterator, string(), ascii::space_type> label;
      qi::rule<Iterator, string(), ascii::space_type> quoted_label;
      qi::rule<Iterator, string(), ascii::space_type> unquoted_label;
      qi::rule<Iterator, double(), ascii::space_type> branch_length;
  };

  template <typename Iterator>
  bool parse_newick(Iterator first, Iterator last, node &root_node)
  {
    newick_grammar<Iterator> grammar;

    bool result = qi::phrase_parse(first, last, grammar, ascii::space, root_node);
    if (first != last) return false;

    return result;
  }

  /***************************** readNewick *********************************/
  /** Reads a tree in Newick format */
  void readNewick (string newick_fn, node &tree) {
    ifstream f_newick;
    f_newick.open(newick_fn.c_str(), ios::in);
    readNewick(f_newick, tree);
    f_newick.close();
  }

  /***************************** readNewick *********************************/
  /** Reads a tree in Newick format */
  void readNewick (istream &input, node &tree) {
    string line;
    getline(input, line);
    if (parse_newick(line.begin(), line.end(), tree) == false) {
      throw runtime_error("Could not parse file in Newick format.");
    }
  }
} // namespace parse
} // namespace treeio