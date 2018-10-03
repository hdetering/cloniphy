#include "evolution.hpp"
#include <cmath>
#include <cstdlib>

using namespace std;

namespace evolution {

/** default c'tor */
GermlineSubstitutionModel::GermlineSubstitutionModel() {
  double q = 1.0/12.0;
  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      if (i==j) {
        Q[i][j] = 0.0;
      }
      else {
        Q[i][j] = q;
      }
    }
  }
}

GermlineSubstitutionModel::GermlineSubstitutionModel (
	const double Qij[4][4]
)
{
  for (int i=0; i<4; ++i)
    for (int j=0; j<4; ++j)
      Q[i][j] = Qij[i][j];
}

GermlineSubstitutionModel::GermlineSubstitutionModel (
	const double p[4], 
	const double titv
)
{
  //kappa=(titv*(p_i[0]+p_i[2])*(p_i[1]+p_i[3]))/(p_i[0]*p_i[2] + p_i[1]*p_i[3]);

  // TODO: add support for GTR model
  // set up HKY mutation rate matrix (A, C, G, T)
  for (int i=0; i<4; ++i)
    Q[i][i] = 0;
  // A -> X
  Q[0][1] = p[1];
  Q[0][2] = p[2]*titv;
  Q[0][3] = p[3];
  // C -> X
  Q[1][0] = p[0];
  Q[1][2] = p[2];
  Q[1][3] = p[3]*titv;
  // G -> X
  Q[2][0] = p[0]*titv;
  Q[2][1] = p[1];
  Q[2][3] = p[3];
  // T -> X
  Q[3][0] = p[0];
  Q[3][1] = p[1]*titv;
  Q[3][2] = p[2];

  // order: A, C, G, T
  /* this does only work with '-std=gnu++11'
  double k = titv;
  Q = {
    {      0,   p[1], k*p[2],   p[3] },
    {   p[0],      0,   p[2], k*p[3] },
    { k*p[0],   p[1],      0,   p[3] },
    {   p[0], k*p[1],   p[2],      0 },
  };*/
}

/** Simulates the nucleotide substitution process for a germline site
  * Assuming that mutation is certain -> use Qij to pick new nucleotide */
short MutateSite(
  short ref_nuc,
  function<double()>& random,
  const GermlineSubstitutionModel& model
)
{
  double r, cumProb[4];
  short new_nuc = 0;

  cumProb[0] = model.Q[ref_nuc][0];
  for (short i=1; i<4; ++i)
    cumProb[i] = cumProb[i-1] + model.Q[ref_nuc][i];
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

/* Nucleotide substition matrix - Jukes-Cantor 1961 */
void GermlineSubstitutionModel::init_JC() {
  double p[4] = { 0.25, 0.25, 0.25, 0.25 };
  double k = 1.0;
  init_HKY(p, k);
}

/* Nucleotide substition matrix - Felseinstein 1981 */
void GermlineSubstitutionModel::init_F81(double p[4]) {
  double k = 1.0;
  init_HKY(p, k);
}

/* Nucleotide substition matrix - Kimura 1980 */
void GermlineSubstitutionModel::init_K80(double k) {
  double p[4] = { 1.0, 1.0, 1.0, 1.0 };
  init_HKY(p, k);
}

/* Nucleotide substition matrix - Hasegawa, Kishino and Yano 1985 */
void GermlineSubstitutionModel::init_HKY(double p[4], double k) {
	kappa = k;
  for (short i=0; i<4; i++)
    for (short j=0; j<4; j++)
      if (i==j) /* nucleotide does not change -> prob 0 (force mutation) */
        Q[i][j] = 0.0;
      else if ((i==0 && j==2) || (i==1 && j==3) || (i==2 && j==0) || (i==3 && j==1)) /* transition */
        Q[i][j] = k*p[j];
      else /* transversion */
        Q[i][j] = p[j];
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


/*----------------------------------------------------------------------------*/
/* SomaticSubstitutionModel
/*----------------------------------------------------------------------------*/

/** Parse TSV file with mutation profile and calculate cumulative frequencies */
SomaticSubstitutionModel::SomaticSubstitutionModel(
    const string& fn_mut_sig,
    const map<string, double>& contrib) :
  m_site(96), m_alt(96), m_weight(96)
{
  // read mutation profiles from files
  int n_lines_read = parseProfiles(fn_mut_sig, contrib);

  // sanity check: did we read complete profiles? expected: 96 probability values
  // (6 types of substitution ∗ 4 types of 5’ base ∗ 4 types of 3’ base)
  if (n_lines_read != 96) {
    fprintf(stderr, "[ERROR] Somatic mutation profiles in file '%s' do not contain the expected number of rows (96)\n", fn_mut_sig.c_str());
  }

  // determine trinuc mutation probs for reverse complement sites
  vector<string> vec_site;
  vector<string> vec_alt;
  vector<double> vec_weight;
  for (auto i=0; i<this->m_site.size(); i++) {
    string site = this->m_site[i];
    string rc_site = seqio::rev_comp(site);
    string alt = this->m_alt[i];
    string rc_alt = seqio::rev_comp(this->m_alt[i]);
    double p = this->m_weight[i];

    vec_site.push_back(rc_site);
    vec_alt.push_back(rc_alt);
    vec_weight.push_back(p);
  }
  this->m_site.insert(this->m_site.end(), vec_site.begin(), vec_site.end());
  this->m_alt.insert(this->m_alt.end(), vec_alt.begin(), vec_alt.end());
  this->m_weight.insert(this->m_weight.end(), vec_weight.begin(), vec_weight.end());
  assert(this->m_site.size() == 192); // sanity check: are all possible mutations represented?
}

/** Read mutation profiles from file and return mutation probabilities for signatures.
 *
 *  Expected input format (tab-separated):
 *  Substitution Type  Trinucleotide  Somatic Mutation Type  Signature 1     Signature 2  ...
 *  C>A                ACA            A[C>A]A                0.011098326166  0.000682708227
 *  C>A                ACC            A[C>A]C                0.009149340734  0.000619107232
 *  ...
 *
 * Input file was be obtained from COSMIC database:
 * http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt
 */
int SomaticSubstitutionModel::parseProfiles(
    const string &filename,
    const map<string, double> &contrib
) {
  int n_lines = 0;
  string line;
  vector<string> header, id_row;
  ifstream filestream;
  filestream.open(filename.c_str());

  // parse header line
  stringio::safeGetline(filestream, line);
  stringio::split(line, '\t', header);

  // read each substitution line, get probabilities for profiles
  stringio::CSVRow row('\t');
  while (filestream >> row && row.size()>0) {
    this->m_site[n_lines] = row[1];
    this->m_alt[n_lines] = row[0].substr(2, 1);
    for (auto i=3; i<row.size(); i++) {
      if (!header[i].empty()) {
        string id_prof = header[i];
        double prob = stod(row[i]);
        if (contrib.find(id_prof) != contrib.end()) {
          this->m_weight[n_lines] += prob * contrib.at(id_prof);
        }
      }
    }

    n_lines++;
  }

  filestream.close();
  return n_lines;
}

} /* namespace evolution */
