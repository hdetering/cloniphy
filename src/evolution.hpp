#ifndef EVOLUTION_H
#define EVOLUTION_H

#include "matrix.hpp"
#include <boost/function.hpp>

/** Provides models of sequence evolution. */
namespace evolution {

struct SubstitutionModel {
  double Pij[4][4]; // per-site substitution rate matrix
  double kappa;   // transition-transversion rate

  // default model: HKY with titv=0.5
  SubstitutionModel(double p_i[4], double titv);
  /** Simulates the nucleotide substitution process for a site */
  short MutateNucleotide(short nuc, boost::function<float()>&);
};

// calculate per-site substitution rates
void HKY (double Pij[4][4], double branchLength, double kappa, double varRate, double p_i[4]);
void GTR (double Pij[4][4], double branchLength, double varRate, double p_i[4]);
void GTnR (double Pij[4][4], double branchLength, double varRate, double p_i[4]);
int EigenREV (double Root[], double Cijk[]);
}

#endif /* EVOLUTION_H */