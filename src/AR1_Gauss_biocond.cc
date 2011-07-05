#include <iostream>
// #include <armadillo>
#include <RcppArmadillo.h>
#include <fstream>
using namespace std;
using namespace arma;
#include "functionsErrorModel.h"
#include "updateFunctions_AR1.h"
#include "commonFunctions.h"
#include <string>
#include <cstdio>
#include <ctime>


void Error_Gauss_c(string& ResultsFolder, mat &x_R, colvec &ParamVec_C)
{
   // Declare variables -----------------------------------------------------------
  // -----------------------------------------------------------------------------
  
  // .. General sampler parameters
  int                                      samples, burnIn, thin;
  double               c, d, a, b, a_exp, b_exp, sigmaS, sigmaMu;
  ucolvec                                      informTimeFlag(2);
  clock_t start_t;
  start_t = clock();
  informTimeFlag.ones();

  // .. Loop vars
  int                        mcmc_iteration;
  
  // .. Strings
  string                 Y_mean_fileName, X_mean_fileName;
  
  // .. Files
  FILE *Bfile, *MuFile, *RhoFile, *LambFile, *GammaFile, *LambExpFile;
  
  // .. MCMC variables and data
  mat                                     Ytplus1, Yt, X, B;
  umat                                             gamma_ij;
  colvec                     lambda_Exp, eta, precTau_i, mu;
  double                                                 ro;
  int                numExp, thinDataSample,  time_m, genes;  
  cube                                              allData;
  mat                                Y_mean,  Xhat, N_it, Y;
  colvec                                           Xhat_sqr;
  
  // .. Aux parameters
  int                  numDiag, p_sqr, free_gammas;
  mat                          Cplus, cplus_aux, C;
  colvec           shape_lamExp, mean_xt, mean_xt1;
  double            logS, eta_mu, eta_s, shape_eta;
  
  // .. Aux variables
  double                   logRosMinlogS, sum_gammas;
  mat      mu_mat, residuals, B_times_xt, precMatrix;
  int                                        counter;
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  
  // .. Read parameter file 
  paramFromVec_Gauss(ParamVec_C, samples, burnIn, thin, c, d, a, b, a_exp, b_exp, 
		     sigmaS, sigmaMu, numExp, thinDataSample);
  
  // .. Open output files
  openOutputFiles_Gauss(ResultsFolder, Bfile, MuFile, RhoFile, LambFile, GammaFile, LambExpFile);
 
  // .. Set random Seed
  Rcpp::RNGScope scope;

  // .. Load data file and scale
  //   int numExp = 4;
  // .. Read data, scale and calulate auxiliary constants. Also init variable Y
  readDataBioreps_ReturnAll(allData, Xhat_sqr, Xhat, N_it, Y, genes, time_m, x_R, numExp);
  
  // .. MCMC variables 
  initMCMCvars_AR1(mu, ro, gamma_ij, B, eta, genes, time_m);  
  // .. No need to init with value, just size
  lambda_Exp.zeros(genes);
  
  counter        = 0;
//   thinDataSample = 6;
  free_gammas    = genes*(genes-1);
  numDiag        = genes;

  // .. Init diagonal elements on
  gamma_ij.diag().ones();
  B = B%gamma_ij;

  // ..  Auxiliary constants
  p_sqr        = genes*genes;
  shape_eta    = a + 0.5 * time_m;
  shape_lamExp = a_exp + 0.5*sum(N_it,1);
  eta_mu       = pow(sigmaMu,-2);
  eta_s        = pow(sigmaS, -2);
  logS         = log(sigmaS);
  precMatrix   = eta_s*eye<mat>(genes, genes);

  Y_mean       = zeros<mat>(genes, time_m+1);
  Yt           = Y.cols( 0,  time_m-1);
  Ytplus1      = Y.cols( 1,    time_m);

  
  mu_mat     = repmat(mu, 1, time_m);
    
  // .. DEBUG
  //   mat B_true;
  //   B_true.load("../DATA/Networks/Lock.TrueNet2");
  //   gamma_ij  = B_true!=0;
  // .. DEBUG
  
  // .. MCMC sampling loop
  for(mcmc_iteration = 0; mcmc_iteration < samples; mcmc_iteration++)
  {  
    if(mcmc_iteration%thinDataSample == 0)
    {
      update_lambdaExp(lambda_Exp, Xhat_sqr, Xhat, N_it, shape_lamExp, b_exp, Y);
      update_Y(Y, Yt, Ytplus1, Xhat, N_it, lambda_Exp, eta, time_m, mu, B, genes);
    }

    // .. Update Rho
    sum_gammas = accu(gamma_ij) - numDiag; 
    ro         = Rf_rbeta(c + sum_gammas, d + free_gammas - sum_gammas);  
    // .. Rho dependant constant
    logRosMinlogS = log(ro/(1-ro)) - logS;

    mean_xt1 = mean(Ytplus1,1);
    mean_xt  = mean(Yt,1);
    // .. Sample from mu
    updateMu_AR1(mu, eta, eta_mu, B, mean_xt1, mean_xt, time_m);

    // .. Aux parameters
    mu_mat = repmat(mu, 1, time_m);
    reCalcYstats(C, Cplus , Yt, Ytplus1, mu_mat);
 
    // .. Sample from eta
    residuals   = B * Yt - Ytplus1 + mu_mat;
    updateEta(eta, residuals, shape_eta, b);    
        
    // .. Update Gibbs variables and regression coefficients
    updateCoeffAndGibbsVars(B, gamma_ij, eta, C, Cplus,  precMatrix, logRosMinlogS, genes);

    if(mcmc_iteration%10 == 0)
      estimateTime_AllowCancel(informTimeFlag, mcmc_iteration, samples, start_t);

    
    if((mcmc_iteration>burnIn) & (mcmc_iteration%thin == 0))
    {
      Y_mean = Y_mean + Y;
      counter++;
      // Timer? Estimated time?
      // Residuals?
      // Write variables to file
      writeMatToFile(Bfile, B);     
      writeToFileDouble(RhoFile, ro);      
      writeToFileVec(LambFile, eta);
      writeToFileVec(LambExpFile, lambda_Exp);
      writeToFileVec(MuFile, mu);
      writeToFileInt(GammaFile, gamma_ij);      
    }    
  }

  fclose(Bfile);
  fclose(RhoFile);
  fclose(LambFile);
  fclose(MuFile);
  fclose(GammaFile);
  fclose(LambExpFile);
  // .. Write inferred posterior data mean
  Y_mean_fileName = ResultsFolder + "Y_mean";  
  Y_mean          = Y_mean/counter;
  Y_mean.save(Y_mean_fileName, raw_ascii);
  X_mean_fileName = ResultsFolder + "DataMean_standarised";  
  Xhat = Xhat/N_it;
  Xhat.save(X_mean_fileName, raw_ascii);
}

