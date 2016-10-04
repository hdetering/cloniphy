#ifndef EVOLUTION_H
#define EVOLUTION_H

#include "matrix.hpp"
#include <functional>

/** Provides models of sequence evolution. */
namespace evolution {

struct SubstitutionModel {
  double Qij[4][4]; // mutation rate matrix (diagonal entries 0)
  double Pij[4][4]; // per-site substitution rate matrix
  double kappa;   // transition-transversion rate

  /** default c'tor */
  SubstitutionModel();
  /**custom transition probabilities */
  SubstitutionModel(double Qij[4][4]);
  // default model: HKY with titv=0.5
  SubstitutionModel(double p_i[4], double titv);
  /** Simulates the nucleotide substitution process for a site */
  short MutateNucleotide(short nuc, std::function<double()>&);
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

} /* namespace evolution */

#endif /* EVOLUTION_H */
