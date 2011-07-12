#include <iostream>
#include <fstream>
#include <iomanip>
// #include <armadillo>
#include <RcppArmadillo.h>
#include <cstdio>
using namespace std;
using namespace arma;

#include <string>
#include "matrixManipulationFunctions.h"
#include "updateFunctions_AR1.h"
#include "commonFunctions.h"


mat nan_cubeMean(const cube& dataMat, const mat& N_it);
cube repMat_cube(const mat& A, int slices);
mat nan_cubeSum(const cube& dataMat);
mat subNaNForZero(mat A);


void openOutputFiles_Gauss(string& ResultsFolder,    FILE* &Bfile,    FILE* &MuFile, FILE* &RhoFile, 
						  FILE* &LambFile, FILE* &GammaFile, FILE* &LambExpFile)
{  
   // .. Declare strings
   string B_name, Mu_name, rho_name, lamb_name, gamma_name, lamb_exp_name; 
   // .. Read numbers into it.
   B_name        = ResultsFolder + "B_mcmc"; 
   Mu_name       = ResultsFolder + "Mu_mcmc"; 
   rho_name      = ResultsFolder + "Rho_mcmc"; 
   lamb_name     = ResultsFolder + "Lambda_mcmc"; 
   gamma_name    = ResultsFolder + "Gamma_mcmc"; 
   lamb_exp_name = ResultsFolder + "Lambda_exp_mcmc"; 

   // .. Open files

   Bfile       = fopen(B_name.c_str(), "w");
   MuFile      = fopen(Mu_name.c_str(), "w");
   RhoFile     = fopen(rho_name.c_str(), "w");
   LambFile    = fopen(lamb_name.c_str(), "w");
   GammaFile   = fopen(gamma_name.c_str(), "w");
   LambExpFile = fopen(lamb_exp_name.c_str(), "w");
}

void paramFromVec_Gauss(const colvec& ParamVec_C, int& samples, int& burnIn, int& thin, double& c, 
		   double& d, double& a, double& b, double& a_exp, double& b_exp, 
		   double& sigmaS, double& sigmaMu, int& numExp, int& thinDataSample)
{  
   samples        = ParamVec_C(0);
   burnIn         = ParamVec_C(1);
   thin           = ParamVec_C(2);
   c              = ParamVec_C(3);
   d              = ParamVec_C(4);
   sigmaS         = ParamVec_C(5);
   a              = ParamVec_C(6);
   b              = ParamVec_C(7);
   a_exp          = ParamVec_C(8);
   b_exp          = ParamVec_C(9);
   sigmaMu        = ParamVec_C(10);   
   thinDataSample = ParamVec_C(11);   
   numExp         = ParamVec_C(12);      
}
void openOutputFiles_Student(string& ResultsFolder,    FILE* &Bfile,    FILE* &MuFile, FILE* &RhoFile, 
						  FILE* &LambFile, FILE* &GammaFile, FILE* &LambExpFile, FILE* &DegFreedomFile)
{  
   // .. Declare strings
   string B_name, Mu_name, rho_name, lamb_name, gamma_name, lamb_exp_name, DegFreedom_name; 
   // .. Read numbers into it.
   B_name          = ResultsFolder + "B_mcmc"; 
   Mu_name         = ResultsFolder + "Mu_mcmc"; 
   rho_name        = ResultsFolder + "Rho_mcmc"; 
   lamb_name       = ResultsFolder + "Lambda_mcmc"; 
   gamma_name      = ResultsFolder + "Gamma_mcmc"; 
   lamb_exp_name   = ResultsFolder + "Lambda_exp_mcmc"; 
   DegFreedom_name = ResultsFolder + "DegFreedom_mcmc"; 

   // .. Open files

   Bfile          = fopen(B_name.c_str(), "w");
   MuFile         = fopen(Mu_name.c_str(), "w");
   RhoFile        = fopen(rho_name.c_str(), "w");
   LambFile       = fopen(lamb_name.c_str(), "w");
   GammaFile      = fopen(gamma_name.c_str(), "w");
   LambExpFile    = fopen(lamb_exp_name.c_str(), "w");
   DegFreedomFile = fopen(DegFreedom_name.c_str(), "w");
}

void paramFromVec_Student(const colvec& ParamVec_C, int& samples, int& burnIn, int& thin, double& c, 
		   double& d, double& a, double& b, double& a_exp, double& b_exp, double& a_deg, 
		   double& b_deg, double& sigmaS, double& sigmaMu, int& numExp, int& thinDataSample)
{  
   samples        = ParamVec_C(0);
   burnIn         = ParamVec_C(1);
   thin           = ParamVec_C(2);
   c              = ParamVec_C(3);
   d              = ParamVec_C(4);
   sigmaS         = ParamVec_C(5);
   a              = ParamVec_C(6);
   b              = ParamVec_C(7);
   a_exp          = ParamVec_C(8);
   b_exp          = ParamVec_C(9);
   a_deg          = ParamVec_C(10);
   b_deg          = ParamVec_C(11);
   sigmaMu        = ParamVec_C(12);   
   thinDataSample = ParamVec_C(13);   
   numExp         = ParamVec_C(14);      
}

// Initialise MCMC variables
void initMCMCvars_Student(colvec &mu, double &ro, umat &gamma_ij, mat &B, colvec &eta, colvec &lambda_Exp, colvec &degFreedom, int genes, int conditions)
{
  double r_min    = 0.0001;
  double r_max    = 0.2;
  double pb_min   = -1;
  double pb_max   = 1;
  double lamb_min = 0.1;
  double lamb_max = 1;
  double student_min = 10;
  double student_max = 20;
  
  // .. Set size of vars
  B.set_size(genes, genes);
  gamma_ij.set_size(genes, genes);
  eta.set_size(genes);
  mu.set_size(genes);
  lambda_Exp.set_size(genes);
  degFreedom.set_size(genes);
  
  // .. Init MCMC variables    
  ro =  Rf_runif( r_min, r_max);  
  RandomBernVec(gamma_ij.memptr(), ro, genes*genes);
  RandomUniformVec(mu.memptr(), pb_min, pb_max, genes);
  RandomUniformVec(B.memptr(), pb_min, pb_max, genes*genes);
  RandomUniformVec(eta.memptr(), lamb_min, lamb_max, genes);
  RandomUniformVec(lambda_Exp.memptr(), student_min, student_max, genes);
  RandomUniformVec(degFreedom.memptr(), student_min, student_max, genes);
}

  
// Read data, subtract mean of each row and divide by stdv
void readDataBioreps_ReturnAll(cube& allData, colvec& Xhat_sqr, mat& Xhat, mat& N_it, mat& Y, int& genes, int& time_m, const mat& dataMat, int numDataSets)
{  
  // .. Declare Vars
  mat                  Xhat_sqr_aux, rep_MeanX, meanMat, stdMat, mean_vec, std_vec; // dataMat, 
  unsigned int                                                          total_time;
  umat                                                                missing_vals;
  
  // .. Load data matrix
//   dataMat.load(filename);
  
  genes      = dataMat.n_rows;
  total_time = dataMat.n_cols/numDataSets;
  time_m     = dataMat.n_cols/numDataSets -1;
  // .. Set data cube size
  allData = cube(dataMat.memptr(), genes, dataMat.n_cols/numDataSets, numDataSets); 
  N_it    = numDataSets*ones<mat>(genes, total_time);
  // .. Calculate N_it number of reps per gene x timepoint
  for(int Dataset_i = 0; Dataset_i < numDataSets; Dataset_i++)
  {
    N_it = N_it - is_NaN_mat(allData.slice(Dataset_i)); // .is_finite();
  }

  rep_MeanX = nan_cubeMean(allData, N_it);
  mean_vec  = mean(rep_MeanX, 1);
  std_vec   = stddev(rep_MeanX, 0, 1);
  
  // .. Subtract mean from all data
  allData = allData - cube(repMat_cube(mat(repmat(mean_vec, 1, total_time)), numDataSets));
  allData = allData / cube(repMat_cube(mat(repmat( std_vec, 1, total_time)), numDataSets));
  
  // .. Calculate Xhat and Xhat_sqr
  Xhat         = nan_cubeSum(allData);
  Xhat_sqr_aux = nan_cubeSum(cube(allData%allData)); //square(allData)
  Xhat_sqr     = sum(Xhat_sqr_aux,1);

  // .. Initialise Y (mcmc variable)
  Y = Xhat/N_it; 
}


void reCalcYstats(mat &C, mat &Cplus , const mat &Yt, const mat &Ytplus1, const mat &mu_mat)
{
  // .. Calculate the Covariance (defined as sum(xt_i*xt_j))
  C     = Yt * trans(Yt);
  // .. and the Cplus matrix (defined as sum(xtplus1_i*xt_j))
  Cplus = (Ytplus1 - mu_mat)*trans(Yt);
}

// ..=========================================
// ..
// ..         Update hidden variables
// ..
// ..=========================================
void update_Y(mat& Y, mat& Yt, mat& Ytplus1, const mat Xhat, const mat N_it, 
	      const colvec& lambda_Exp, const colvec& eta, int time_m, const colvec& mu, 
	      const mat& B, int genes)
{
  mat         B_eta_B, inv_covMat, covMatY;
  colvec                    mu_Y, aux_mu_Y;
  int                                    t;
  vec                         diag_B_eta_B;
  
  // .. Aux vars 
  B_eta_B      = trans(B)*diagmat(eta)*B;
  diag_B_eta_B = B_eta_B.diag();
  
  // .. Calculate distribution parameters
  inv_covMat = B_eta_B + diagmat(N_it.col(0)%lambda_Exp);
  covMatY    = inv(inv_covMat);
  aux_mu_Y   = trans(trans((Y.col(1) - mu)%eta) * B) + Xhat.col(0)%lambda_Exp;
  mu_Y       = covMatY*aux_mu_Y;
  symmetriseMat(covMatY);  
//     covMatY    = (covMatY + trans(covMatY))/2;
  
  // .. Sample
  Y.col(0)    = mvnrnd(mu_Y, covMatY);

  for(t = 1; t < time_m; t++)
  {
    // .. Calculate distribution parameters
    inv_covMat.diag() = diag_B_eta_B + eta + N_it.col(t)%lambda_Exp;
//         inv_covMat = B_eta_B + diagmat(eta + N_it.col(t)%lambda_Exp);
    covMatY    = inv(inv_covMat);
    aux_mu_Y   = trans(trans((Y.col(t+1) - mu)%eta) * B) + (mu + B*Y.col(t-1))%eta + Xhat.col(t)%lambda_Exp;
    mu_Y       = covMatY*aux_mu_Y;
    symmetriseMat(covMatY);  
//     covMatY    = (covMatY + trans(covMatY))/2;
    // .. Sample  
    Y.col(t)   = mvnrnd(mu_Y, covMatY);
  }
  
  t = time_m;
  // .. Calculate distribution parameters
  inv_covMat = diagmat(eta + N_it.col(t)%lambda_Exp);
  covMatY    = inv(inv_covMat);
  aux_mu_Y   = (mu + B*Y.col(t-1))%eta + Xhat.col(t)%lambda_Exp;
  mu_Y       = covMatY*aux_mu_Y;
    symmetriseMat(covMatY);  
//     covMatY    = (covMatY + trans(covMatY))/2;
  // .. Sample  
  Y.col(time_m) = mvnrnd(mu_Y, covMatY);

  Yt      = Y.cols(0, time_m -1);
  Ytplus1 = Y.cols(1, time_m);
}

void update_lambdaExp(colvec &lambda_Exp, const colvec &Xhat_sqr, 
		      const mat &Xhat, const mat &N_it, const colvec &shape_lamExp, double b_exp, const mat& Y)
{
  colvec YsqrxN, sum_squares, rate_eta, scale_eta;
  
  YsqrxN          = sum(mat(square(Y)%N_it),1);
  sum_squares     = YsqrxN + Xhat_sqr - 2*sum(mat(Xhat%Y),1);
  rate_eta        = b_exp + 0.5 * sum_squares;
  scale_eta       = 1./rate_eta;
  // .. Sample
  for(unsigned int loop_var = 0; loop_var < rate_eta.n_elem ;loop_var++)
  {
    lambda_Exp[loop_var] = Rf_rgamma( shape_lamExp[loop_var], scale_eta[loop_var]);
  }
}

mat nan_cubeMean(const cube& dataMat, const mat& N_it)
{
  mat outMean;
  outMean = nan_cubeSum(dataMat);
/*  outMean.zeros(N_it.n_rows, N_it.n_cols);
  for(unsigned int Dataset_i = 0; Dataset_i < dataMat.n_slices; Dataset_i++)
  {
    outMean = outMean + dataMat.slice(Dataset_i);
  }*/
  outMean = outMean/N_it;
  return(outMean);
}

cube repMat_cube(const mat& A, int slices)
{
  // .. Declare and size cube
  cube cube_out(A.n_rows, A.n_cols, slices);
  // Fill
  for(int Dataset_i = 0; Dataset_i < slices; Dataset_i++)
    cube_out.slice(Dataset_i) = A;
  // Return
  return(cube_out);
}

mat nan_cubeSum(const cube& dataMat)
{
  mat outSum;
  outSum.zeros(dataMat.n_rows, dataMat.n_cols);
  for(unsigned int Dataset_i = 0; Dataset_i < dataMat.n_slices; Dataset_i++)
  {
    outSum = outSum + mat(subNaNForZero(dataMat.slice(Dataset_i)));
  }
  return(outSum);
}

mat subNaNForZero(mat A)
{
  mat::iterator A_iter;
  mat::iterator first_elem = A.begin();
  mat::iterator  last_elem = A.end();
  for( A_iter = first_elem; A_iter != last_elem; A_iter++)
  {
    if(isnan(*A_iter))
    {
	*A_iter = 0;
    }
  } 
  return(A);
}

void update_MH_DegFreedom_t(colvec &degFreedom, ucolvec &numAccepted, double a_deg, double b_deg, const cube &w_itr, int genes, 
			    const colvec &numSamples_i, double a_RW)
{
  colvec                                       b_RW_inv, newDegVals(genes), b_RW_inv_prime;
  colvec                                          alpha, subtract_degs, half_old, half_new; 
  colvec           logRatioProposal(genes), logGammaDeg_new(genes), logGammaDeg_old(genes);
  colvec     logHastings_Ratio, logRatio_aux1, logRatio_aux2, logRatio_aux3, logRatio_post;
  colvec                                                               randUnif_vec(genes);
  mat                                                                            alpha_mat;
  ucolvec isAccepted;  
  
  // .. Gamma Random walk
  b_RW_inv  = degFreedom/a_RW;
  
  // .. Sample new proposed values
  for(unsigned int loop_var = 0; loop_var < b_RW_inv.n_elem; loop_var++)
  {
    newDegVals[loop_var] = Rf_rgamma(a_RW, b_RW_inv[loop_var]);
  }

  b_RW_inv_prime = newDegVals/a_RW;
  // .. Calculate log ratio of the proposal
  for(int loop_var = 0; loop_var < genes; loop_var++)
  {
    // Use log pdf option of dgamma
    logRatioProposal[loop_var] = Rf_dgamma(degFreedom[loop_var], a_RW, b_RW_inv_prime[loop_var], 1) - Rf_dgamma(newDegVals[loop_var], a_RW, b_RW_inv[loop_var], 1);
  }

  // .. Aux variables 
  subtract_degs   = newDegVals - degFreedom;
  half_old        = 0.5*degFreedom;
  half_new        = 0.5*newDegVals;
  for(unsigned int loop_var = 0; loop_var < b_RW_inv.n_elem; loop_var++)
  {
    logGammaDeg_new[loop_var] = Rf_lgammafn(half_new[loop_var]);
    logGammaDeg_old[loop_var] = Rf_lgammafn(half_old[loop_var]);
  }

  // Posterior log ratio
  logRatio_aux1 = (a_deg - 1)*(log(newDegVals)-log(degFreedom));
  logRatio_aux2 = subtract_degs%( 0.5 * sum(mat(nan_cubeSum(cube(log(w_itr)))), 1) - 0.5*sum(mat(nan_cubeSum(w_itr)),1) - b_deg);
  logRatio_aux3 = numSamples_i%(half_new%colvec(log(half_new))-half_old%colvec(log(half_old))-logGammaDeg_new+logGammaDeg_old);
  logRatio_post = logRatio_aux1 + logRatio_aux2 + logRatio_aux3;

  // Hastings log ratio
  logHastings_Ratio = logRatio_post + logRatioProposal;
  alpha_mat         = join_rows(zeros<colvec>(genes), logHastings_Ratio);
  alpha             = min(alpha_mat, 1);
  
  // .. Which have been accepted
  RandomUniformVec(randUnif_vec.memptr(), 0., 1., genes);  
  isAccepted  = alpha > colvec(log(randUnif_vec));
  numAccepted = numAccepted + isAccepted;
  // .. Change appropriate values
  placeVecInVec_FromVec(degFreedom, newDegVals, isAccepted);
}


void update_weights_t(cube &w_itr, cube &YminX_sqr, const mat &Y, const cube &dataSet, 
		      const colvec &degFreedom, const mat &lambda_exp, int reps, int time_m)
{
  cube                        Y_rep, half_degs_rep, half_lambda_exp_rep;
  cube                               shape_weights, rate_eta, scale_eta;

  Y_rep               = repMat_cube(Y, reps);
  YminX_sqr           = square(cube(Y_rep - dataSet));
  half_degs_rep       = repMat_cube(mat(repmat(degFreedom/2., 1, time_m + 1)), reps);
  half_lambda_exp_rep = repMat_cube(mat(repmat(lambda_exp/2., 1, time_m + 1)), reps);
  shape_weights       = half_degs_rep + 0.5;
  rate_eta            = half_degs_rep + half_lambda_exp_rep % YminX_sqr;
  scale_eta           = 1./rate_eta;
  
  // .. Sample
  for(unsigned int loop_var = 0; loop_var < rate_eta.n_elem; loop_var++)
  {
    w_itr[loop_var] = Rf_rgamma(shape_weights[loop_var], scale_eta[loop_var]);
  }

}

void update_Y_tDist(mat& Y, mat& Yt, mat& Ytplus1, const cube& dataSets, const cube& w_itr, 
		    const colvec& lambda_Exp, const colvec& eta, int time_m, const colvec& mu, 
		    const mat& B, int genes)
{
  mat         B_eta_B, inv_covMat, covMatY, weightedSum, sumWeights;
  colvec                                             mu_Y, aux_mu_Y;
  int                                                             t;
  vec                                                  diag_B_eta_B;
   
  // .. Aux vars 
  B_eta_B      = trans(B)*diagmat(eta)*B;
  weightedSum  = nan_cubeSum(cube(w_itr%dataSets));
  sumWeights   = nan_cubeSum(w_itr);
  diag_B_eta_B = B_eta_B.diag();

  // .. Calculate distribution parameters
  inv_covMat = B_eta_B + diagmat(sumWeights.col(0)%lambda_Exp);
  covMatY    = inv(inv_covMat);
  aux_mu_Y   = trans(trans((Y.col(1) - mu)%eta) * B) + weightedSum.col(0)%lambda_Exp;
  mu_Y       = covMatY*aux_mu_Y;
  symmetriseMat(covMatY);  
//     covMatY    = (covMatY + trans(covMatY))/2;
  
  // .. Sample
  Y.col(0)    = mvnrnd(mu_Y, covMatY);

  for(t = 1; t < time_m; t++)
  {
    // .. Calculate distribution parameters
    inv_covMat.diag() = diag_B_eta_B + eta + sumWeights.col(t)%lambda_Exp;

//     inv_covMat = B_eta_B + diagmat(eta + sumWeights.col(t)%lambda_Exp);
    covMatY    = inv(inv_covMat);
    aux_mu_Y   = trans(trans((Y.col(t+1) - mu)%eta) * B) + (mu + B*Y.col(t-1))%eta + weightedSum.col(t)%lambda_Exp;
    mu_Y       = covMatY*aux_mu_Y;
    symmetriseMat(covMatY);  
//     covMatY    = (covMatY + trans(covMatY))/2;
    // .. Sample  
    Y.col(t)   = mvnrnd(mu_Y, covMatY);
  }
    
  t = time_m;
  // .. Calculate distribution parameters
  inv_covMat = diagmat(eta + sumWeights.col(t)%lambda_Exp);
  covMatY    = inv(inv_covMat);
  aux_mu_Y   = (mu + B*Y.col(t-1))%eta + weightedSum.col(t)%lambda_Exp;
  mu_Y       = covMatY*aux_mu_Y;
  symmetriseMat(covMatY);  
//     covMatY    = (covMatY + trans(covMatY))/2;
  // .. Sample  
  Y.col(time_m) = mvnrnd(mu_Y, covMatY);

  Yt      = Y.cols(0, time_m -1);
  Ytplus1 = Y.cols(1, time_m);
}


void update_LambdaExp_t(colvec &lambda_exp, const cube &YminX_sqr,  const cube &w_itr, const colvec &shape_lamExp, double b_exp)
{
  colvec rate_eta, scale_lambda_exp;
  
  rate_eta         = b_exp +  0.5*sum(nan_cubeSum(cube(w_itr % YminX_sqr)),1);
  scale_lambda_exp = 1./rate_eta;
  // .. Sample
  for(unsigned int loop_var = 0; loop_var < rate_eta.n_elem; loop_var++)
  {
    lambda_exp[loop_var] = Rf_rgamma(shape_lamExp[loop_var], scale_lambda_exp[loop_var]);
  }
}
