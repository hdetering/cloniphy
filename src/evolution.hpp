#ifndef EVOLUTION_H
#define EVOLUTION_H

#include "matrix.hpp"

/** Provides models of sequence evolution. */
namespace evolution {

struct SubstitutionModel {
  double Pij;
  double kappa;
  // default model: HKY with titv=0.5
  SubstitutionModel(double p_i[4], double titv=0.5);

  void HKY (double Pij[4][4], double branchLength, double kappa, double varRate, double p_i[4]);
  void GTR (double Pij[4][4], double branchLength, double varRate, double p_i[4]);
  void GTnR (double Pij[4][4], double branchLength, double varRate, double p_i[4]);
};

int EigenREV (double Root[], double Cijk[]);
}

#endif /* EVOLUTION_H */
