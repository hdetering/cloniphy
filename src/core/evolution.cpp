#include "evolution.hpp"
#include <cmath>
#include <cstdlib>

namespace evolution {

/** default c'tor */
SubstitutionModel::SubstitutionModel() {
  double q = 1.0/12.0;
  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      if (i==j) {
        Qij[i][j] = 0.0;
      }
      else {
        Qij[i][j] = q;
      }
    }
  }
}

/**custom transition probabilities */
SubstitutionModel::SubstitutionModel(double Q[4][4]) {
  for (int i=0; i<4; ++i)
    for (int j=0; j<4; ++j)
      this->Qij[i][j] = Q[i][j];
}

SubstitutionModel::SubstitutionModel(double p[4], double titv) {
  //kappa=(titv*(p_i[0]+p_i[2])*(p_i[1]+p_i[3]))/(p_i[0]*p_i[2] + p_i[1]*p_i[3]);

  // TODO: add support for GTR model
  // set up HKY mutation rate matrix (A, C, G, T)
  for (int i=0; i<4; ++i)
    Qij[i][i] = 0;
  // A -> X
  Qij[0][1] = p[1];
  Qij[0][2] = p[2]*titv;
  Qij[0][3] = p[3];
  // C -> X
  Qij[1][0] = p[0];
  Qij[1][2] = p[2];
  Qij[1][3] = p[3]*titv;
  // G -> X
  Qij[2][0] = p[0]*titv;
  Qij[2][1] = p[1];
  Qij[2][3] = p[3];
  // T -> X
  Qij[3][0] = p[0];
  Qij[3][1] = p[1]*titv;
  Qij[3][2] = p[2];

  // order: A, C, G, T
  /* this does only work with '-std=gnu++11'
  double k = titv;
  Qij = {
    {      0,   p[1], k*p[2],   p[3] },
    {   p[0],      0,   p[2], k*p[3] },
    { k*p[0],   p[1],      0,   p[3] },
    {   p[0], k*p[1],   p[2],      0 },
  };*/
}

/** Simulates the nucleotide substitution process for a site
  * Assuming that mutation is certain -> use Qij to pick new nucleotide */
short SubstitutionModel::MutateNucleotide(short ref_nuc, std::function<double()>& random) {
	double r, cumProb[4];
  short new_nuc = 0;

  cumProb[0] = Qij[ref_nuc][0];
	for (short i=1; i<4; ++i)
		cumProb[i] = cumProb[i-1] + Qij[ref_nuc][i];
  // normalize
  for (short i=0; i<4; ++i)
    cumProb[i] /= cumProb[3];

	r = random();
	if (r >= 0.0 && r <= cumProb[0])
		new_nuc = 0;
	else if (r > cumProb[0] && r <= cumProb[1])
		new_nuc = 1;
	else if (r > cumProb[1] && r <= cumProb[2])
		new_nuc = 2;
	else
		new_nuc = 3;

  return new_nuc;
}


/* parameters */
double pinv;      // proportion of invariable sites
double titv;      // transition/transversion ratio
double Rmat[6];   // transition rate matrix A-C A-G A-T C-G C-T G-T, for GTR models
double NRmat[12]; // general rate matrix (AC CA AG GA AT TA CG GC CT TC GT=1 TG)

/* global variables */
double Pij;        // per-site substitution rate matrix
double Qij[16];    // rate matrix
double Cijk[256];
double Root[4];
double mr;
double tstv;

// set default values for parameters
void init() {
  pinv = 0.0; // proportion of invariable sites
  // transition rate matrix A-C A-G A-T C-G C-T G-T, for GTR models
  for (int i = 0; i < 6; i++)	Rmat[i] = -1;
  // general rate matrix (AC CA AG GA AT TA CG GC CT TC GT=1 TG)
  for (int i = 0; i < 12; i++) NRmat[i] = -1;
}


/*---------------------------------- HKY -------------------------------------*/
/**	HKY performs Hasegawa-Kishino-Yano 85 correction */
void HKY (double Pij[4][4], double branchLength, double kappa, double varRate, double p_i[4])
{
	int			i, j;
	double	A, t, PIj, beta;

	beta = 0.5 / ((p_i[0] + p_i[2])*(p_i[1] + p_i[3]) + kappa*((p_i[0]*p_i[2]) + (p_i[1]*p_i[3])));

	if (varRate > 0)
		varRate = varRate / (1.0 - pinv);

	t = branchLength * varRate;
	/*fprintf(stderr,"\n t = %lf \n", t);*/

	if (t == 0.0)  	/* no mutations */
		for (i=0; i<4; i++)	{
			for (j=0; j<4; j++)	{
				if (i == j)
					Pij[i][j] = 1.0;
				else
					Pij[i][j] = 0.0;
			}
		}
	else				/* there are mutations */
		for (i=0; i<4; i++)	{
		  for (j=0; j<4; j++) {
				if (j == 0 || j == 2)	/* purine */
					PIj = p_i[0] + p_i[2];
				else
					PIj = p_i[1] + p_i[3]; /* pyrimidine */

				A = 1 + PIj*(kappa-1);

				if (i==j) /* diagonal principal */
					Pij[i][j] = p_i[j] + p_i[j]*(1/PIj - 1)*exp(-beta*t) + ((PIj-p_i[j])/PIj)*exp(-beta*t*A);
				else if ((i==0 && j==2) || (i==1 && j==3) || (i==2 && j==0) || (i==3 && j==1)) /* transition */
					Pij[i][j] = p_i[j] + p_i[j]*(1/PIj - 1)*exp(-beta*t) - (p_i[j]/PIj)*exp(-beta*t*A);
				else /* transversion */
					Pij[i][j] = p_i[j]*(1-exp(-beta*t));
			}
		}
}

/*---------------------------------- GTR -------------------------------------*/
/** General Time-Reversable model */
void GTR (double Pij[4][4], double branchLength, double varRate, double p_i[4])
{
	int 	i, j, k;
	double	t, expt[4];
	/* double Rmat[6], Qij[16], Cijk[256], Root[4], mr, tstv; Global Variables*/

	if (varRate > 0)
		varRate = varRate / (1.0 - pinv);

	t = branchLength * varRate;

	k=0;
	for (i=0; i<3; i++)
		for (j=i+1; j<4; j++)
      if (i*4+j != 11)
				Qij[i*4+j]=Qij[j*4+i]=Rmat[k++];

	Qij[3*4+2]=Qij[2*4+3]=1.0;

	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			Qij[i*4+j] *= p_i[j];

	mr=0;
	for (i=0; i<4; i++) {
		Qij[i*4+i]=0;
		Qij[i*4+i]=-(Qij[i*4]+Qij[i*4+1]+Qij[i*4+2]+Qij[i*4+3]);

		mr-=p_i[i]*Qij[i*4+i];
	}

	EigenREV(Root, Cijk);

	/* calculate mean ts/tv ratio */
	mr=2*(p_i[3]*Qij[3*4+1]+p_i[0]*Qij[0*4+2]);
	tstv=mr/(1-mr);

	/* P(t)ij = SUM Cijk * exp{Root*t}*/
	if (t<1e-6)
		{
		for (i=0; i<4; i++)
			for (j=0; j<4; j++)
				{
				if (i==j)
					Pij[i][j] = 1.0;
				else
					Pij[i][j] = 0.0;
				}
		}
	else
		{
		for (k=1; k<4; k++)
			expt[k]=exp(t*Root[k]);
		for (i=0; i<4; i++)
			for (j=0; j<4; j++)
				{
				Pij[i][j]=Cijk[i*4*4+j*4+0];
				for (k=1; k<4; k++)
					Pij[i][j]+=Cijk[i*4*4+j*4+k]*expt[k];
				}
		}
}

/*----------------------------------- GTnR -----------------------------------*/
/* GTR non reversible */
void GTnR (double Pij[4][4], double branchLength, double varRate, double p_i[4])
{
	int 	i, j, k;
	double	t, expt[4];

	if (varRate > 0)
		varRate = varRate / (1.0 - pinv);

	t = branchLength * varRate;

	/*
		A	C	G	T
	A	0	1	2	3
	C	4	5	6	7
	G	8	9	10	11
	T	12	13	14	15
	*/

	/* fills no symmetrical matrix */
	k=0;
	for (i=0; i<3; i++)
		for (j=i+1; j<4; j++) {
			Qij[i*4+j]=NRmat[k++];
			Qij[j*4+i]=NRmat[k++];
		}

/*	AC CA AG GA AT TA CG GC CT TC GT=1 TG */

	/* all rates relative to GT */
	Qij[11] = 1.0;

	for (i=0; i<4; i++)
		for (j=0; j<4; j++)
			Qij[i*4+j] *= p_i[j];

	mr=0;
	for (i=0; i<4; i++) {
		Qij[i*4+i]=0;
		Qij[i*4+i]=-(Qij[i*4]+Qij[i*4+1]+Qij[i*4+2]+Qij[i*4+3]);
		mr-=p_i[i]*Qij[i*4+i];
	}

	EigenREV(Root, Cijk);

/* calculate mean ts/tv ratio */ /*double check*/
/*	mr=2*(p_i[3]*Qij[3*4+1]+p_i[0]*Qij[0*4+2]);*/
	mr = p_i[3]*Qij[3*4+1] + p_i[0]*Qij[0*4+2] + p_i[1]*Qij[1*4+3] + p_i[2]*Qij[2*4+0] ;
	tstv=mr/(1-mr);

/* P(t)ij = SUM Cijk * exp{Root*t}
*/
	if (t<1e-6) { /* too small branch */
		for (i=0; i<4; i++)
			for (j=0; j<4; j++)	{
				if (i==j)
					Pij[i][j] = 1.0;
				else
					Pij[i][j] = 0.0;
			}
	}
	else {
		for (k=1; k<4; k++)
			expt[k]=exp(t*Root[k]);
		for (i=0; i<4; i++)
			for (j=0; j<4; j++) {
				Pij[i][j]=Cijk[i*4*4+j*4+0];
				for (k=1; k<4; k++)
					Pij[i][j]+=Cijk[i*4*4+j*4+k]*expt[k];
			}
	}
}

/* Eigen function for nucleotide models */
int EigenREV (double Root[], double Cijk[])
{
	int i,j,k;
	double U[16], V[16], T1[16], T2[16];

	matrix::abyx (1/mr, Qij, 16);

	if ((k=matrix::eigen (1, Qij, 4, Root, T1, U, V, T2))!=0) {
		fprintf(stderr, "\ncomplex roots in EigenREV");
		exit(0);
	}
	matrix::xtoy (U, V, 16);
	matrix::matinv (V, 4, 4, T1);
	for (i=0; i<4; i++)
   		for (j=0; j<4; j++)
   			for (k=0; k<4; k++)
   				Cijk[i*4*4+j*4+k] = U[i*4+k]*V[k*4+j];

	return (0);
}

} /* namespace evolution */
