#include <iostream>
#include <fstream>
#include <iomanip>
// #include <armadillo>
#include <RcppArmadillo.h>
#include <cstdio>

using namespace std;
using namespace arma;
#include <string>
#include "commonFunctions.h"
#include "matrixManipulationFunctions.h"



void initBasesOn(urowvec& bases_on, const umat& gamma_ij, int i, int M);
mat bspline_mat(const colvec& x, const double&xl, const double&xr, const int& ndx, const int& deg);
void MHStep_Splines( urowvec& bases_on,   urowvec& gamma_Row,  double& logMVPDF_Old, int i, int j, 
	    const mat& lambxCPlusS,   const rowvec& lambxCplusIdot, const colvec& sumLogs, int M);
void calc_logMVPDF_withBases(double& logMVPDF, const mat& lambxCPlusS , const rowvec& lambxCplusIdot, const unsigned int& i, urowvec& gamma_Row);
void updateCoefficients_splines(                mat& B,              const int& i, const urowvec& links_on, 
			const mat& lambxCPlusS, const rowvec& lambxCplusIdot);
double rTruncGamma(double oldVal, double a_shape, double b_scale, double truncate_me);
void modifyBasesOnVector(urowvec& vec, int j, int M, int newVal);


void openOutputFiles_Splines(string& ResultsFolder, FILE* &MuFile, FILE* &LambdaFile, FILE* &DesignFile, FILE* &GammaFile, 
		     FILE* &RhoFile, FILE* &TauFile, FILE* &All_fFile)
{  
   // .. Declare strings
   string B_name, Mu_name, rho_name, lamb_name, gamma_name, design_name, tau_name, allF_name; 
   // .. Read numbers into it.
   Mu_name      = ResultsFolder + "Mu_mcmc"; 
   lamb_name    = ResultsFolder + "Lambda_mcmc"; 
   design_name  = ResultsFolder + "SplineDesignMat"; 
   gamma_name   = ResultsFolder + "Gamma_mcmc";    
   rho_name     = ResultsFolder + "Rho_mcmc"; 
   tau_name     = ResultsFolder + "Tau_mcmc"; 
   allF_name     = ResultsFolder + "MeanPosterior_Functions"; 
   
   // .. Open files
   MuFile     = fopen(Mu_name.c_str(), "w");
   LambdaFile = fopen(lamb_name.c_str(), "w");
   RhoFile    = fopen(rho_name.c_str(), "w");
   GammaFile  = fopen(gamma_name.c_str(), "w");
   TauFile    = fopen(tau_name.c_str(), "w");
//    DesignFile = fopen(design_name.c_str(), "w");
//    All_fFile  = fopen(allF_name.c_str(), "w");

}

void paramFromVec_Splines(const colvec& ParamVec_C, int& samples, int& burnIn, int& thin, double& c, 
		   double& d, double& a, double& b, double& truncate_me, double& tau0, int& M,
		   double& sigmaMu, double& a_pareto)
{  
   samples     = ParamVec_C(0);
   burnIn      = ParamVec_C(1);
   thin        = ParamVec_C(2);
   c           = ParamVec_C(3);
   d           = ParamVec_C(4);
   truncate_me = ParamVec_C(5);
   tau0        = ParamVec_C(6);
   M           = ParamVec_C(7);
   a           = ParamVec_C(8);
   b           = ParamVec_C(9);   
   sigmaMu     = ParamVec_C(10);
   a_pareto    = ParamVec_C(11);
}


// .. Build the B-spline design mat for all regressions
mat despline(const mat& data, const int& ndx, const int& deg)
{
  // .. Declare vars
  colvec                   x;
  int             g, M, T, G;
  mat          bs, designMat;
  double              xl, xr;
  
  // .. Genes and timepoints
  G = data.n_cols;
  T = data.n_rows; 
  M = ndx + deg;

  designMat.zeros(T, G*M);
  for(g = 0; g < G; g++)
  {
      x  = data.col(g);
      xl = min(x) - 0.1; 
      xr = max(x) + 0.1;
      bs = bspline_mat(x,xl,xr,ndx,deg);
      designMat.cols(g*M, (g+1)*M - 1) = bs;    
  //     sum(bs,2) 
  }
  return(designMat);
}


// .. Calculate design matrix for single regression
mat bspline_mat(const colvec& x, const double&xl, const double&xr, const int& ndx, const int& deg)
{
  double             dx;
  mat        X, P, B, B_aux, T;
  rowvec     aux_seq, t;
  colvec              r;
  umat               aux_logic_mat3,   aux_logic_mat1, aux_logic_mat2;

  dx = (xr-xl)/ndx;
  aux_seq = generate_seq(-deg, ndx-1);
  t  = xl + dx *aux_seq;
 
  T  = repmat(t, x.n_elem, 1); 
  X  = repmat(x, 1, t.n_elem);  
  P  = (X-T)/dx;
  aux_logic_mat1 = (T<=X);
  aux_logic_mat2 = (X < (T+dx));
  aux_logic_mat3 = umat(aux_logic_mat1 + aux_logic_mat2) == 2;
  B = conv_to<mat>::from(aux_logic_mat3);
  
  // .. Build index Order Vector
  r.set_size(t.n_elem);
  r.rows(0, t.n_elem - 2) = trans(generate_seq(1, t.n_elem - 1));
  r(r.n_elem - 1) = 0;
  
  for(int k = 1; k < (deg+1); k++)
  {
    B_aux = reorderMatColsFromVec(B, r);
    B = (P%B + (k+1-P)%B_aux)/k;
  }
  return(B);
}

// .. Build matrix that maps gamma index to b index
umat buildMapGammaBeta(int M, int genes)
{
  rowvec            aux_seq;
  umat         mapGammaBeta;
  aux_seq = generate_seq(0, M*genes-1);
  aux_seq.reshape(M, genes);
  mapGammaBeta = conv_to<umat>::from(aux_seq);
  return(mapGammaBeta);
}

// .. Init B with Gamma
void fillBzerosUseGamma(mat& B, const umat& gamma_ij, int M)
{
  urowvec bases_on(gamma_ij.n_cols*M);
  for(unsigned int i=0; i<gamma_ij.n_cols; i++)
  {
    initBasesOn(bases_on, gamma_ij, i, M);
    B.row(i) = bases_on % B.row(i);
  }
}

void initBasesOn(urowvec& bases_on, const umat& gamma_ij, int i, int M)
{
  // .. Vars
  urowvec                  links_on;
  int       G, start_indx, end_indx;
  
  // .. Vector to loop through
  links_on = gamma_ij.row(i);
  G        = links_on.n_elem;
  for(int indx = 0; indx < G; indx++)
  {
    start_indx = M*indx;
    end_indx   = M*(indx+1)-1;
    if(links_on[indx])
    {
      modifyBasesOnVector(bases_on, indx, M, 1);
//       bases_on.cols(start_indx, end_indx).fill(1);
    }else{
      modifyBasesOnVector(bases_on, indx, M, 0);
//       bases_on.cols(start_indx, end_indx).fill(0);      
    }
  }
}

// .. Use number of bases M ans index j to place newVal in vec
void modifyBasesOnVector(urowvec& vec, int j, int M, int newVal)
{
  int start_loop = j*M;
  int end_loop   = (j+1)*M;
  for(int loop_var = start_loop; loop_var < end_loop; loop_var++ )
  {
    vec[loop_var] = newVal;
  }
}

// Initialise MCMC variables
void initMCMCvars_Splines(colvec &mu, double &ro, umat &gamma_ij, mat &B, colvec &eta, int genes, int conditions, mat &tau_ij, int M)
{
  double r_min    = 0.0001;
  double r_max    = 0.2;
  double pb_min   = -1;
  double pb_max   = 1;
  double lamb_min = 0.1;
  double lamb_max = 1;
  
  // .. Set size of vars
  B.set_size(genes, M*genes);
  tau_ij.set_size(genes, genes);
  gamma_ij.set_size(genes, genes);
  eta.set_size(genes);
  mu.set_size(genes);
  
  // .. Init MCMC variables    
  ro =  Rf_runif(r_min, r_max);  
  RandomBernVec(gamma_ij.memptr(), ro, genes*genes);
  RandomUniformVec(mu.memptr(), pb_min, pb_max, genes);
  RandomUniformVec(B.memptr(), pb_min, pb_max, genes*genes*M);
  RandomUniformVec(eta.memptr(), lamb_min, lamb_max, genes);
  RandomUniformVec(tau_ij.memptr(), lamb_min, lamb_max, genes*genes);

}


// .. Build precision of B prior 
mat priprec(int M)
{
  mat P(M,M);
  P.zeros();
  
  for(int j=2; j < (M-2); j++)
  {
      P(j,j)=6;
      P(j,j-1)= -4;   P(j,j+1)= -4;
      P(j-1,j)= -4;   P(j+1,j)= -4;
      P(j,j-2) = 1;   P(j,j+2) = 1;
      P(j-2,j) = 1;   P(j+2,j) = 1;
  }
  
  P(M-1,M-1)=1;
  P(M-1,M-2)=-2; P(M-2,M-1)=-2;
  P(M-1,M-3)=1;  P(M-3,M-1)=1;

  P(M-2,M-2)=5;
  P(M-2,M-3)=-4;  P(M-3,M-2)=-4;
  P(M-2,M-4)=1;  P(M-4,M-2)=1;

  P(0,0)= 1; 
  P(0,1)=-2; P(1,0)=-2;
  P(1,2)=-4;  P(2,1)=-4;

  P(1,3) = 1; P(3,1)=1;
  P(1,1) = 5; 
  
  return(P);
}

void updateTaus(mat& tau_ij, colvec& logRosMinlogS, const mat& smallPriorMat, const umat& gamma_ij, const mat& B, 
		double shape_tau_self, double shape_tau, int M, double logRoDivOneMinRo, double truncate_me, double a_pareto_inv,
		double m_minus2_div2, int i, double tau0)
{
  // .. Vars 
  colvec                  logS(tau_ij.n_cols);
  double      uni_rnd, rate_tau, scale_tau, b_Prior_prec_b, shape_tau_me;
  urowvec                            links_on;
  rowvec                            b_ij, b_i;
  int                               num_genes;
  num_genes =  (int) tau_ij.n_cols;
  links_on  = gamma_ij.row(i);
  b_i       = B.row(i);
  // .. With i fixed loop over j and sample from Tau_ij
  for(int j_tau = 0; j_tau < num_genes; j_tau++)
  {
    // .. Set shape parameter to one of two (self or not)
    if(i == j_tau)
    {
      shape_tau_me = shape_tau_self;
    }else{
      shape_tau_me = shape_tau;
    }
    // .. If connection is on: sample as usual. If not sample from prior
    if(links_on[j_tau])
    {
      b_ij             = b_i.cols(M * j_tau, M  * (j_tau + 1) - 1);
      b_Prior_prec_b   = as_scalar(b_ij * smallPriorMat * trans(b_ij));
      rate_tau         = 0.5 * b_Prior_prec_b;
      scale_tau        = 1/rate_tau;
      tau_ij(i, j_tau) = rTruncGamma(tau_ij(i, j_tau), shape_tau_me, scale_tau, truncate_me);
    }else{
      uni_rnd = unif_rand();
      tau_ij(i, j_tau) = truncate_me*pow(uni_rnd, a_pareto_inv);
    }
    // DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//     tau_ij(i,j_tau) = 500;
    logS(j_tau) = -log(tau0 * pow(tau_ij(i,j_tau), m_minus2_div2));
  }      
  logRosMinlogS  = logRoDivOneMinRo - logS;       
}

// .. Sample from a truncated gamma distribution
double rTruncGamma(double oldVal, double a_shape, double b_scale, double truncate_me)
{
  double      uni_rnd, x_new, y, a_inv, k;
  bool               flag_BadSample;
  // .. To speed things up select between methods using simple criterion (mode vs truncation point)
  if((a_shape-1)*b_scale < truncate_me)
  {
    // .. Rejection sampling =====================================================
    flag_BadSample = true;
    while(flag_BadSample)
    {
      // .. Sample from gamma
      x_new  = Rf_rgamma(a_shape, b_scale);
      // .. Check sampled value comes from truncated region
      if (x_new < truncate_me)
      {
	flag_BadSample   = false;
      } 
    }
  }else{
    // .. Damien and Walker 2001 (Very fast, but heavy autocorrelation) ==========
    y       =  Rf_runif(0, exp(-oldVal/b_scale));
    a_inv   = 1/a_shape;
    k       = min_two(truncate_me, -log(y)*b_scale);
    uni_rnd = unif_rand();
    x_new   = pow(uni_rnd,a_inv)*k;
  }
  return(x_new);
}

// .. Build prior matrix
void priorMultiTau(mat& priorWithTau, const mat& smallPriorMat, const rowvec& tau_i, double tau0, int M, int genes)
{  
  // .. Vars
  int   start_indx, end_indx;  
  mat                  subMe;
  // .. Fill mat
  for(int block_j = 0; block_j < genes;  block_j++) 
  { 
    // .. Calc big mat indices 
    start_indx = block_j*M;
    end_indx   = (block_j + 1)*M - 1;
    // .. Multiply plain mat times tau
    subMe      = smallPriorMat*tau_i[block_j];
    // .. Incorporate tau0
    subMe(0,0) = subMe(0,0) + tau0;
    subMe(1,1) = subMe(1,1) + tau0;
    // .. Substitute
    priorWithTau.submat(start_indx, start_indx, end_indx, end_indx) = subMe; 
  }    
}

// ..   Sample from  Gammas      
void updateGammaAndB_row_i(mat& B,   umat& gamma_ij, const mat& FullprecB, const rowvec& meanBNoPrec, 
			     const colvec & logRosMinlogS, int genes, int M, int i)
{
  // .. Vars
  urowvec     bases_on(M*genes), links_on;
  double                     logMVPDF_Old;
  int                                   j;
  ucolvec                   seq_rnd(genes);

  
  // .. Links in row i that are on and corresponding bases
  links_on    = gamma_ij.row(i);
  initBasesOn(bases_on, gamma_ij, i, M);
  calc_logMVPDF_withBases(logMVPDF_Old, FullprecB, meanBNoPrec, i, bases_on);      

  // .. Random permutation of update index
  random_intSequence(seq_rnd);
  // .. Loop over elements of row i and sample from gamma
  for(int j_loop = 0; j_loop < genes; j_loop++)
  {
    // .. Get permuted index
    j = seq_rnd[j_loop];
    // .. MH update link i,j We reuse logMVPDF_Old (value 
    // .. for current state of gamma_row)
    // .. Self interactions are not updated
    if(i!= j)
    {
      MHStep_Splines(bases_on, links_on, logMVPDF_Old, i, j, FullprecB, meanBNoPrec, logRosMinlogS, M);
    }
  }      
  // .. Update gamma matrix with new row
  gamma_ij.row(i)= links_on; 
  
  // .. Update B
  updateCoefficients_splines(B, i, bases_on, FullprecB, meanBNoPrec);
}



// .. Sample from mu  -----------------------------------------------
void updateMu_Splines(colvec& mu, const colvec& eta, const double& eta_mu, const mat& B, const colvec& mean_xt1, const rowvec& xhat_1, int time_m, int i)
{
    // Declare vars
    double aux_denominator, aux_muVar, mu_sd, mu_mean;
     
    // .. Aux variables for mu
    aux_denominator = eta(i)*time_m;    
    aux_muVar = 1/(1 + eta_mu / aux_denominator);    
    mu_sd     = sqrt(aux_muVar/aux_denominator);
    mu_mean   = aux_muVar * mean_xt1(i) - xhat_1(i);       
        
    // .. Sample
    mu(i) = Rf_rnorm(mu_mean, mu_sd);
}

// .. C and X are NOT passed by reference
void updateCoefficients_splines(                mat& B,              const int& i, const urowvec& links_on, 
			const mat& lambxCPlusS, const rowvec& lambxCplusIdot)
{
  // .. Declare vars
  mat                     covarianceMatrixB,  lambxCPlusSReduced;
  rowvec        lambxCplusIdotReduced, lambxCplusIdotReduced_aux;
  unsigned int                                            num_on;
  colvec                                         aux_mu_B, b_new;
  
  // .. Before updating links, check if any need updating
  num_on   = accu(links_on);
  if (num_on>0)
  {
      // .. Aux variables for B    
      subMatFromVector(lambxCPlusSReduced, lambxCPlusS, links_on);
//       lambxCplusIdotReduced_aux = lambxCplusIdot.row(i);
      subVectorFromVector(lambxCplusIdotReduced, lambxCplusIdot, links_on);
      // .. B Covariance matrix
      covarianceMatrixB = inv(lambxCPlusSReduced);
      
      // .. Correct precision => symmetry problem
      covarianceMatrixB = (covarianceMatrixB + trans(covarianceMatrixB))/2;
      // .. Update B for links that are ON
      aux_mu_B  = covarianceMatrixB*trans(lambxCplusIdotReduced);
      b_new     = mvnrnd(aux_mu_B, covarianceMatrixB);
  }
  // .. Update B row with new values and zeros
  fillMatRowWithVecAndZeros(B, b_new, i, links_on);
}

void calcIndivF(mat& all_f, mat& all_f_sqr, mat& full_F_sqr, const mat& Y, const mat& B, int genes, int M, int time_m)
{
  colvec             B_i, spline_function_ij;
  mat           indiv_f(time_m, genes*genes);
  int                   start_indx, end_indx;
  
  // .. In order to write out individual splines
  for(int i = 0; i < genes; i++)
  {  
    B_i = trans(B.row(i));
    for(int j = 0; j < genes; j++)
    {
      start_indx               = j*M;
      end_indx                 = (j+1)*M - 1;     
      spline_function_ij       = Y.cols(start_indx, end_indx) * B_i.rows(start_indx, end_indx);
      indiv_f.col(genes*i + j) = spline_function_ij - mean(spline_function_ij); 
    }
    full_F_sqr.col(i) = full_F_sqr.col(i) + square(sum(indiv_f.cols(genes*i, genes*(i+1) - 1), 1));
  }
  all_f     = all_f + indiv_f;
  all_f_sqr = all_f_sqr + square(indiv_f);
}



// .. Metropolis-Hastings move to update a gamma_ij
// .. NOTE: Using log porbabilities. Always calculate move from
// .. link off to on. If move is inverese then multiply by -1
void MHStep_Splines( urowvec& bases_on,   urowvec& gamma_Row,  double& logMVPDF_Old, int i, int j, 
	    const mat& lambxCPlusS,   const rowvec& lambxCplusIdot, const colvec& sumLogs, int M)
{
  // .. Vars
  unsigned int               gamma_old, gamma_proposed;  
  double         signOfK, alfa, log_aux, hastingsRatio; 
  double          logMVPDF_1, logMVPDF_0, logMVPDF_new;
  colvec                         bugVec;// , hastingsRatio
  colvec zero_colvec(1);
  zero_colvec(0) =0;
  // .. Get current value of gamma_ij  
  gamma_old    = gamma_Row(j);
  // .. Change current value of gamma_ij (new proposed config)
  if (gamma_old){ 
    gamma_proposed = 0; 
  }else{ 
    gamma_proposed = 1;
  }
  gamma_Row(j) = gamma_proposed;
  
//   bases_on.cols(j*M, (j+1)*M-1).fill(gamma_proposed);
  modifyBasesOnVector(bases_on, j, M, gamma_proposed);

  
  // .. Calculate logMVPDF for this link config
  calc_logMVPDF_withBases(logMVPDF_new, lambxCPlusS, lambxCplusIdot, i, bases_on);      
  
   // .. Which is logMVPDF for gamma_ij=1 and which for gamma_ij=0
  if (gamma_old){    
      logMVPDF_1 = logMVPDF_Old;
      logMVPDF_0 = logMVPDF_new;
      signOfK    = -1;
  }
  else{
      logMVPDF_1 = logMVPDF_new;
      logMVPDF_0 = logMVPDF_Old;
      signOfK    = 1;
  }

  // .. Calculate hastings ratio for move from gamma_ij
  // .. and change sign if signOfK is -1
  hastingsRatio = sumLogs[j] + 0.5 *(logMVPDF_1 - logMVPDF_0);
  hastingsRatio = signOfK * hastingsRatio;

    alfa          = min_two(0., hastingsRatio);
/*  hastingsRatio.insert_rows(0,zero_colvec);
  alfa          = min(hastingsRatio);*/
  log_aux       = log(unif_rand());
/*  if (i == 1)
  {
    cout << hastingsRatio << endl;  
    cout << alfa << endl << "---------------" << endl;  
  }*/
  // .. To accept new move
  if (alfa > log_aux){
      // .. Move has been accepted, no need to modify gamma_Row.
      // .. Old value of logMVPDF (for next iteration)
      logMVPDF_Old = logMVPDF_new;  
  }
  else
  {
      // .. Move has been rejected change gamma_ij back to old       
      // .. value logMVPDF_Old is still the same
      gamma_Row(j) = gamma_old;
      modifyBasesOnVector(bases_on, j, M, gamma_old);
//       bases_on.cols(j*M, (j+1)*M-1).fill(gamma_old);
  } 
}


void calc_logMVPDF_withBases(double& logMVPDF, const mat& lambxCPlusS , const rowvec& lambxCplusIdot, const unsigned int& i, urowvec& gamma_Row)
{
  // .. Vars
  unsigned int                                               num_on;  
  mat                                            lambxCPlusSReduced;
  rowvec           lambxCplusIdotReduced_aux, lambxCplusIdotReduced;
  
 
  // .. Make sure there is at least one link
  num_on   = accu(gamma_Row);
  if(num_on>0){	
      // .. Calculate logMVPDF for current state of gamma_row
      subMatFromVector(lambxCPlusSReduced, lambxCPlusS, gamma_Row);
      subVectorFromVector(lambxCplusIdotReduced, lambxCplusIdot, gamma_Row);
      MHlogMVPDF(logMVPDF, lambxCPlusSReduced, lambxCplusIdotReduced);
  }
  else{
      logMVPDF = 0;
  }
 
}


