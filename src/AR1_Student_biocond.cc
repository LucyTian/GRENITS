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
// ZZ Profiler!!!!!
// #include <google/profiler.h>

void Error_Student_c(string& ResultsFolder, mat &x_R, colvec &ParamVec_C, mat &Gamma_fixed)
{
  // Declare variables -----------------------------------------------------------
  // -----------------------------------------------------------------------------
  
  // .. General sampler parameters
  int                                            samples, burnIn, thin;
  double       a_deg, b_deg, c, d, a, b, a_exp, b_exp, sigmaS, sigmaMu;
  ucolvec                                            informTimeFlag(2);
  clock_t start_t;
  informTimeFlag.ones();

  // .. Loop vars
  int                        mcmc_iteration;
  
  // .. Strings
  string                                 param_filName, Y_mean_fileName;
  string                          X_mean_fileName, NumAccepted_fileName;
  
  // .. Files
  FILE *Bfile, *MuFile, *RhoFile, *LambFile, *GammaFile, *LambExpFile, *DegFreedomFile;
  
  // .. MCMC variables and data
  mat                                     Ytplus1, Yt, X, B;
  umat                                             gamma_ij;
  colvec         degFreedom, lambda_Exp, eta, precTau_i, mu;
  double                                                 ro;
  int                numExp, thinDataSample,  time_m, genes;  
  cube                            YminX_sqr, w_itr, allData;
  mat                                Y_mean,  Xhat, N_it, Y;
  colvec                                           Xhat_sqr;
  
  // .. Aux parameters
  int           num_fixedON, numDiag, p_sqr, free_gammas;
  mat                                           Cplus, C;
  colvec   numSamples_i, shape_lamExp, mean_xt, mean_xt1;
  double                  logS, eta_mu, eta_s, shape_eta;
  double                                     CV_RW, a_RW;

  
  // .. Aux variables
  double                       logRosMinlogS, sum_gammas;
  mat          mu_mat, residuals, B_times_xt, precMatrix;
  int                            counter_accept, counter;
  ucolvec                                    numAccepted;
  umat                                  regMat, UpdateMe;
  ucolvec   flatRegsIndx_Vec, UpdateIndx_Vec, numRegsVec;
   
  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  // ZZ Profiler!!!!!
//   ProfilerStart("../../ProfileStudent");

  // .. Read parameter file 
  paramFromVec_Student(ParamVec_C, samples, burnIn, thin, c, d, a, b, a_exp, b_exp, 
		       a_deg, b_deg, sigmaS, sigmaMu, numExp, thinDataSample);
  
  // .. Open output files
  openOutputFiles_Student(ResultsFolder, Bfile, MuFile, RhoFile, LambFile, GammaFile, LambExpFile, DegFreedomFile);
   
  // .. Set random Seed
  Rcpp::RNGScope scope;

  // .. Read data, scale and calulate auxiliary constants. Also init variable Y
  readDataBioreps_ReturnAll(allData, Xhat_sqr, Xhat, N_it, Y, genes, time_m, x_R, numExp);

  // .. Save mean data set
  X_mean_fileName = ResultsFolder + "DataMean_standarised";  
  Xhat = Xhat/N_it;
  Xhat.save(X_mean_fileName, raw_ascii);

  // .. Free memory of unused variables
  Xhat_sqr.reset();
  Xhat.reset();
  
  // .. MCMC variables 
  initMCMCvars_Student(mu, ro, gamma_ij, B, eta, lambda_Exp, degFreedom, genes, time_m);  

  // .. Process Gamma_fixed matrix: Fix gammas, calc parameters
  processFixedGammas(Gamma_fixed, num_fixedON,    free_gammas, UpdateMe, gamma_ij, 
		      numRegsVec,      regMat, UpdateIndx_Vec, flatRegsIndx_Vec);

  w_itr.zeros(genes, time_m+1, numExp);  
  counter        = 0;
  counter_accept = 0;
  numAccepted    = zeros<ucolvec>(genes);
  numSamples_i   = sum(N_it,1);


  // .. Tuning parameter for gamma random walk sampler
  CV_RW = 0.3; 
  a_RW  = pow(CV_RW,-2);

  
  // .. Init diagonal elements on
  gamma_ij.zeros();
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
  
  start_t = clock();

  // .. MCMC sampling loop
  for(mcmc_iteration = 1; mcmc_iteration < samples+1; mcmc_iteration++)
  {  
    if(mcmc_iteration%thinDataSample == 0)
    {
      update_weights_t(w_itr, YminX_sqr, Y, allData, degFreedom, lambda_Exp, numExp, time_m);
      update_LambdaExp_t(lambda_Exp, YminX_sqr,  w_itr, shape_lamExp, b_exp);
      update_MH_DegFreedom_t(degFreedom, numAccepted, a_deg, b_deg, w_itr, genes, numSamples_i, a_RW);   
//       degFreedom << 13.266923 << 9.311586 << 4.491057 << 21.695926 << 5.877218 << endr;
      update_Y_tDist(Y, Yt, Ytplus1, allData, w_itr, lambda_Exp, eta, time_m, mu, B, genes);      
      counter_accept++;
    }
//     update_MH_DegFreedom_t(degFreedom, numAccepted, a_deg, b_deg, w_itr, genes, numSamples_i, a_RW);   

    // .. Update Rho
    sum_gammas = accu(gamma_ij) - num_fixedON; 
    ro         = Rf_rbeta( c + sum_gammas, d + free_gammas - sum_gammas);  
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
    updateCoeffAndGibbsVars_reg(            B, gamma_ij,      eta,          C,  Cplus,  eta_s, 
				logRosMinlogS,    genes, UpdateMe, numRegsVec, regMat);
//     updateCoeffAndGibbsVars(B, gamma_ij, eta, C, Cplus,  precMatrix, logRosMinlogS, genes);

    if(mcmc_iteration%10 == 0)
      estimateTime_AllowCancel(informTimeFlag, mcmc_iteration, samples, start_t);

    
    if((mcmc_iteration>burnIn) & (mcmc_iteration%thin == 0))
    {
      Y_mean = Y_mean + Y;
      counter++;
      // Timer? Estimated time?
      // Residuals?
      // Write variables to file
      writeMatToFile_withIndx(Bfile, B, flatRegsIndx_Vec);     
      writeToFileDouble(RhoFile, ro);      
      writeToFileVec(LambFile, eta);
      writeToFileVec(LambExpFile, lambda_Exp);
      writeToFileVec(DegFreedomFile, degFreedom);      
      writeToFileVec(MuFile, mu);
      writeToFileInt_withIndx(GammaFile, gamma_ij, UpdateIndx_Vec);      
    }    
  }
  fclose(Bfile);
  fclose(RhoFile);
  fclose(LambFile);
  fclose(MuFile);
  fclose(GammaFile);
  fclose(LambExpFile);
  fclose(DegFreedomFile);
  
  // .. Write inferred posterior data mean
  Y_mean_fileName = ResultsFolder + "Y_mean";  
  Y_mean          = Y_mean/counter;
  Y_mean.save(Y_mean_fileName, raw_ascii);
  NumAccepted_fileName = ResultsFolder + "acceptanceRatio";  
  colvec dub_numAccepted = conv_to<colvec>::from(numAccepted)/counter_accept;
  dub_numAccepted.save(NumAccepted_fileName, raw_ascii);  
  
  // ZZ Profiler!!!!!
//   ProfilerStop();

}

