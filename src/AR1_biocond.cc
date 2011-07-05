#include <iostream>
// #include <armadillo>
#include <RcppArmadillo.h>

#include <fstream>
using namespace std;
using namespace arma;
#include "updateFunctions_AR1.h"
#include "commonFunctions.h"

#include <string>
#include <cstdio>
#include <ctime>


void AR1_c(string& ResultsFolder, mat &x_R, colvec &ParamVec_C)
{
  // Declare variables -----------------------------------------------------------
  // -----------------------------------------------------------------------------
  
  // .. General sampler parameters
  int                        samples, burnIn, thin; //, debugInt;
  double               c, d, a, b, sigmaS, sigmaMu;
  ucolvec                        informTimeFlag(2);
  clock_t start_t;
  start_t = clock();
  informTimeFlag.ones();

  // .. Loop vars
  int                        mcmc_iteration;
  
  // .. Strings
  string  dataFile, param_filName; //ResultsFolder,
  
  // .. Files
  FILE *Bfile, *MuFile, *RhoFile, *LambFile, *GammaFile;
  
  // .. MCMC variables and data
  mat                     xt_plus1, xt, X, B;
  umat                              gamma_ij;
  colvec                  eta, precTau_i, mu;
  double                                  ro;
  int          time_m, genes, num_Conditions;  
  
  // .. Aux parameters
  int            numDiag, p_sqr, free_gammas;
  mat                    Cplus, cplus_aux, C;
  colvec                   mean_xt1, mean_xt;
  double      logS, eta_mu, eta_s, shape_eta;
  
  // .. Aux variables
  double                   logRosMinlogS, sum_gammas;
  mat      mu_mat, residuals, B_times_xt, precMatrix;
  
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  
  // .. Read parameter file 
  paramFromVec_AR1(ParamVec_C, samples, burnIn, thin, c, d, a, b, sigmaS, sigmaMu);
  
  // .. Open output files
  openOutputFiles_AR1(ResultsFolder, Bfile, MuFile, RhoFile, LambFile, GammaFile);

  X = ScaleData(x_R);
    
  //   X.load(dataFile.c_str());
  genes          = X.n_rows;      
  num_Conditions = X.n_cols;
  time_m         = num_Conditions - 1;
  
  // Destroy X ? zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
  xt       = X.cols( 0,  time_m-1);
  xt_plus1 = X.cols( 1,    time_m);

  // .. Set random Seed  
  Rcpp::RNGScope scope;

  
  // .. MCMC variables 
  initMCMCvars_AR1(mu, ro, gamma_ij, B, eta, genes, num_Conditions);
  
  // .. Init diagonal elements on
  gamma_ij.diag().ones();
  B = B%gamma_ij;

//   Load true coeff values   
//  B.load("../../DATA/MNI_All2Perturbed.TrueNetwork");
  free_gammas = genes*(genes-1);
  numDiag     = genes;

  // .. Auxiliary constants
  p_sqr       = genes*genes;
  shape_eta   = a + 0.5 * time_m;
  eta_mu      = pow(sigmaMu,-2);
  eta_s       = pow(sigmaS, -2);
  mean_xt1    = mean(xt_plus1,1);
  mean_xt     = mean(xt,1);
  logS        = log(sigmaS);
  precMatrix  = eta_s*eye<mat>(genes, genes);
  mu_mat     = repmat(mu, 1, time_m);
  
  //  .. Calculate the Covariance (defined as sum(xt_i*xt_j))
  C         = xt * trans(xt);
  cplus_aux = xt_plus1 * trans(xt);

//   cout << C;
  // .. DEBUG
  /*  mat B_true;
  B_true.load("../DATA/Networks/Lock.TrueNet2");
  gamma_ij  = B_true!=0;*/
  // .. DEBUG
  
  // .. MCMC sampling loop
  for(mcmc_iteration = 0; mcmc_iteration < samples; mcmc_iteration++)
  {  
    // .. Update Rho
    sum_gammas = accu(gamma_ij) - numDiag; 
    ro         = Rf_rbeta(c + sum_gammas, d + free_gammas - sum_gammas);  
    
    // .. Rho dependant constant
    logRosMinlogS = log(ro/(1-ro)) - logS;

    // .. Sample from eta
    residuals   = xt_plus1 - B * xt - mu_mat;   
    updateEta(eta, residuals, shape_eta, b);    
    
    // .. Sample from mu
    updateMu_AR1(mu, eta, eta_mu, B, mean_xt1, mean_xt, time_m);

    // .. Aux parameters
    mu_mat = repmat(mu, 1, time_m);
    Cplus  = cplus_aux - mu_mat*trans(xt);
    
    // .. Update Gibbs variables and regression coefficients
    updateCoeffAndGibbsVars(B, gamma_ij, eta, C, Cplus,  precMatrix, logRosMinlogS, genes);

    if(mcmc_iteration%10 == 0)
      estimateTime_AllowCancel(informTimeFlag, mcmc_iteration, samples, start_t);

    if((mcmc_iteration>burnIn) & (mcmc_iteration%thin == 0))
    {
      // Residuals?
      // Write variables to file
      writeMatToFile(Bfile, B);     
      writeToFileDouble(RhoFile, ro);      
      writeToFileVec(LambFile, eta);
      writeToFileVec(MuFile, mu);
      writeToFileInt(GammaFile, gamma_ij);      
    }    
  }    
  fclose(Bfile);
  fclose(RhoFile);
  fclose(LambFile);
  fclose(MuFile);
  fclose(GammaFile);
/*  return 0;*/
}

