#ifndef DATAFRAME_H
#define DATAFRAME_H

#include "../stringio.hpp"

//#include <cstdlib> /* atof */
#include <cassert>
#include <string>
#include <vector>

/** Collection of useful data classes. */
namespace model {

/** Convert data element to string. */
template <typename T>
std::string data_to_str(const T data);

/** Convert string to data element. */
template <typename T>
T str_to_data(const std::string);

/** A data matrix with named rows and columns. */
template <typename T>
struct DataFrame {
  /** Number of rows. */
  unsigned n_rows;
  /** Number of columns. */
  unsigned n_cols;
  /** Stores column names. */
  std::vector<std::string> colnames;
  /** Stores row names. */
  std::vector<std::string> rownames;
  /** Stores actual data. */
  std::vector<std::vector<T>> data;

  /** Default c'tor. */
  //DataFrame () : n_rows(0), n_cols(0) {}
  DataFrame ();

  /** C'tor to create DataFrame from matrix of strings.
   *  Assumes that for column to contain row labels and 
   *  first row to contain column labels.
   */
  DataFrame(const std::vector<std::vector<std::string>> data);

  /** Access data elements by label. */
  T at (
    const std::string lbl_row, 
    const std::string lbl_col
  ) const
  {
    unsigned i = 0;
    unsigned j = 0;

    // search row names for label
    auto it_row = this->rownames.begin();
    while (it_row != this-> rownames.end() && *it_row != lbl_row) {
      i++;
      it_row++;
    }
    // search col names for label
    auto it_col = this->colnames.begin();
    while (it_col != this->colnames.end() && *it_col != lbl_col) {
      j++;
      it_col++;
    }
    // make sure labels are valid
    assert ( it_row != this->rownames.end() );
    assert ( it_col != this->colnames.end() );

    return this->data[i][j];
  }

  /** Access data elements by index. */
  T at (
    const unsigned idx_row, 
    const unsigned idx_col
  ) const
  {
    // make sure labels are valid
    assert ( idx_row < this->n_rows );
    assert ( idx_col < this->n_cols );

    return this->data[idx_row][idx_col];
  }

  /** Format as string. */
  std::string to_string() const
  {
    std::string str;
    // print header row
    if (colnames.size() > 0) {
      for (unsigned j=0; j<colnames.size(); j++) {
        str += stringio::format("\t%s", colnames[j].c_str());
      }
      str += "\n";
    }
    // print data rows
    for (unsigned i=0; i<n_rows; i++) {
      str += rownames[i];
      for (unsigned j=0; j<n_cols; j++) {
        str += data_to_str<T>(this->at(i,j));
      }
      str += "\n";
    }

    return str;
  }
};

/* Method implementations */

template <>
inline std::string data_to_str<double>(const double value) {
  return stringio::format("\t%.2f", value);
}

template <>
inline double str_to_data<double>(const std::string str) {
  return stod(str);
}

template <typename T>
DataFrame<T>::DataFrame () 
: n_rows(0), n_cols(0) 
{}

template <typename T>
DataFrame<T>::DataFrame (
  const std::vector<std::vector<std::string>> d
) 
{
  std::vector<std::string> lbl_row;
  std::vector<std::string> lbl_col;
  std::vector<std::vector<T>> mtx_data;

  this->n_rows = d.size() - 1;
  this->n_cols = d[0].size() - 1;

  // read column labels
  for (unsigned j=1; j<n_cols+1; j++) {
    lbl_col.push_back(d[0][j]);
  }
  // read data rows
  for (unsigned i=1; i<n_rows+1; i++) {
    // add row label
    lbl_row.push_back(d[i][0]);
    std::vector<T> data_row;
    // read elements of current row
    for (unsigned j=1; j<n_cols+1; j++) {
      data_row.push_back( str_to_data<T>(d[i][j]) );
    }
    // add new data row
    mtx_data.push_back( data_row );
  }

  this->rownames = lbl_row;
  this->colnames = lbl_col;
  this->data = mtx_data;
}

} /* namespace model */

#endif /* DATAFRAME_H */