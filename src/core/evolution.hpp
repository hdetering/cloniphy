#ifndef EVOLUTION_H
#define EVOLUTION_H

#include "matrix.hpp"
#include "stringio.hpp"
#include <functional>

/** Provides models of germline sequence evolution. */
namespace evolution {

struct GermlineSubstitutionModel {
  double Qij[4][4]; // mutation rate matrix (diagonal entries 0)
  double Pij[4][4]; // per-site substitution rate matrix
  double kappa;   // transition-transversion rate

  /** default c'tor */
  GermlineSubstitutionModel();
  /** custom transition probabilities */
  GermlineSubstitutionModel(double Qij[4][4]);
  /** default model: HKY with titv=0.5 */
  GermlineSubstitutionModel(double p_i[4], double titv);
  // calculate nucleotide substitution matrix (Q/mu)
  void init_JC();
  void init_F81(double p[4]);
  void init_K80(double k);
  void init_HKY(double p[4], double k);
};

// calculate per-site substitution rates
void HKY (double Pij[4][4], double branchLength, double kappa, double varRate, double p_i[4]);
void GTR (double Pij[4][4], double branchLength, double varRate, double p_i[4]);
void GTnR (double Pij[4][4], double branchLength, double varRate, double p_i[4]);
int EigenREV (double Root[], double Cijk[]);

/** Provides access to somatic mutation profiles. */
struct SomaticSubstitutionModel {
  std::vector<std::string> m_site; // list of tri-nucleotides
  std::vector<std::string> m_alt; // list of alternative alleles
  std::vector<double> m_weight; // probabilities for all possible substitutions

  /** default c'tor */
  SomaticSubstitutionModel () {};
  /** Create substitution matrix from COSMIC file. */
  SomaticSubstitutionModel (
      const std::string &filename,
      const std::map<std::string, double> &contrib
  );
  /** Read mutation profiles from file and return cumulative mutation probabilities for signatures. */
  int parseProfiles (
      const std::string &filename,
      const std::map<std::string, double> &contrib);
};

/** Simulates the nucleotide substitution process for a germline site */
short MutateSite(
  short idx_nuc,
  std::function<double()>&,
  const GermlineSubstitutionModel&
);


} /* namespace evolution */

#endif /* EVOLUTION_H */
