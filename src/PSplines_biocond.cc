#include <iostream>
// #include <armadillo>
#include <RcppArmadillo.h>
#include <fstream>
using namespace std;
using namespace arma;
#include "updateFunctions_Splines.h"
#include "matrixManipulationFunctions.h"
#include "commonFunctions.h"
#include <string>
#include <cstdio>
#include <ctime>

void PSplines_c(string& ResultsFolder, mat &x_R, colvec &ParamVec_C)
{
  // Declare variables -----------------------------------------------------------
  // -----------------------------------------------------------------------------
  
  // .. General sampler parameters
  int                                   M, samples, burnIn, thin; //, debugInt;
  double        c, d, a, b, truncate_me, tau0, sigmaMu, a_pareto;
  ucolvec                                      informTimeFlag(2);
  clock_t start_t;
  start_t = clock();
  informTimeFlag.ones();

  // .. Loop vars
  int                        mcmc_iteration;
  
  // .. Strings
  string       full_F_sqr_file, dataFile, param_filName, allf_file, allf_sqr_file;
  
  // .. Files
  FILE *MuFile, *LambdaFile, *DesignFile, *GammaFile, *RhoFile, *TauFile, *All_fFile;
  
  // .. MCMC variables and data
  mat             tau_ij, xt_plus1, xt, X, B;
  umat                              gamma_ij;
  colvec                  eta, precTau_i, mu;
  double                                  ro;
  int          time_m, genes, num_Conditions;  
  
  // .. Aux parameters
  int                nodesSpline, degreeSpline, numDiag, p_sqr, free_gammas;
  double                                                      m_minus2_div2;
  mat                         smallPriorMat, xhat_0, Y, Cplus, cplus_aux, C;
  colvec                                                  mean_xt1, mean_xt;
  rowvec                                                             xhat_1;
  double         shape_tau, shape_tau_self, a_pareto_inv, eta_mu, shape_eta;

  // .. Aux variables
  double                                   logRoDivOneMinRo, sum_gammas;
  mat                YB_full, mu_mat, residuals, B_times_xt, precMatrix;
  mat       full_F_sqr ,priorMatB, FullprecB, all_f, indiv_f, all_f_sqr;
  ucolvec                                                      bases_on;
  colvec                                        YB, logRosMinlogS, logS;
  rowvec                                             tau_i, meanBNoPrec;
  umat                                                     mapGammaBeta;
  int                                                           counter;

  // -----------------------------------------------------------------------------
  // -----------------------------------------------------------------------------
  
  // .. Read parameter file 
  paramFromVec_Splines(ParamVec_C, samples, burnIn, thin, c, d, a, b, 
 		truncate_me, tau0, M, sigmaMu, a_pareto);

  // .. Open output files
  openOutputFiles_Splines(ResultsFolder, MuFile, LambdaFile, DesignFile, GammaFile, RhoFile, TauFile, All_fFile);
  
  // .. Load data file and scale
  X = ScaleData(x_R);

  //   X.load(dataFile.c_str());
  genes          =           X.n_rows;      
  num_Conditions =           X.n_cols;
  time_m         = num_Conditions - 1;
  
  // Destroy X ? zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
  xt       = X.cols( 0,  time_m-1);
  xt_plus1 = X.cols( 1,    time_m);

  // .. Set random Seed
  Rcpp::RNGScope scope;
  
  // .. Build design matrix Y
  degreeSpline = 3;
  nodesSpline  = M - degreeSpline;
  Y = despline(trans(xt), nodesSpline, degreeSpline);
  
  // .. Write Y to file
//   Y.save("DesignFile", raw_ascii);
  
  // .. MCMC variables 
  initMCMCvars_Splines(mu, ro, gamma_ij, B, eta, genes, num_Conditions, tau_ij, M);

  // .. Init diagonal elements on
  gamma_ij.diag().ones();
  
  // .. Init aux variables
  bases_on = zeros<ucolvec>(genes*M);
  
  // .. Create mapGammaBeta matrix, mapping gamma index to B indeces
  fillBzerosUseGamma(B, gamma_ij, M);
   
  // .. Auxiliary constants
  p_sqr          = genes*genes;
  shape_eta      = a + 0.5 * time_m;
  eta_mu         = pow(sigmaMu,-2);
  mean_xt1       = mean(xt_plus1,1);
  m_minus2_div2  = (M-2)/2.0;
  a_pareto_inv   = 1.0/a_pareto;
  shape_tau      = m_minus2_div2 + a_pareto;
  shape_tau_self = m_minus2_div2 + 3.*a_pareto;
  free_gammas    = genes*(genes-1);
  numDiag        = genes;
  
  // .. x hats
  YB_full = Y*trans(B);
  xhat_1  = mean(YB_full);
  xhat_0  = YB_full - repmat(xhat_1, time_m, 1);
  // .. Remove mat
  YB_full.reset();
  
  // .. Calculate the "Covariance"
  C = trans(Y) * Y;
  
  // .. Calculate prior precision matrix (Same for all regressions)
  // .. Prior precision auxiliary constant matrices
  smallPriorMat  = priprec(M);
  priorMatB.zeros(M*genes, M*genes);
//   noTau_PriorMat = DiagnalBlockMat(smallPriorMat, genes);
  
  // .. Init stuff for mean spline function output 
  counter = 0;
  full_F_sqr.zeros(time_m, genes);
  all_f.zeros(time_m, genes*genes);
  all_f_sqr.zeros(time_m, genes*genes);
  // DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/*  gamma_ij.load("../../../DATA/Networks/Lock.TrueNet2");
  gamma_ij.diag().ones();*/
//   gamma_ij.ones();
  // DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // .. MCMC sampling loop
  for(mcmc_iteration = 0; mcmc_iteration < samples; mcmc_iteration++)
  {  

    // .. Update Rho
    sum_gammas = accu(gamma_ij) - numDiag; 
    ro         = Rf_rbeta( c + sum_gammas, d + free_gammas - sum_gammas);  
      // DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//     ro =0.45;
    // .. Rho dependant constant
    logRoDivOneMinRo = log(ro/(1-ro));

    // .. Sample from eta
    mu_mat = repmat(mu, 1, time_m);
    residuals = xt_plus1 - trans(xhat_0) - mu_mat;      
    updateEta(eta, residuals, shape_eta, b); 
    // DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/*    
    eta.fill(1000);
    mu.zeros();
*/
    // .. Loop through genes
    for(int i = 0; i < genes; i++)
    {
      // .. Sample from taus
      updateTaus(tau_ij, logRosMinlogS, smallPriorMat, gamma_ij, B, shape_tau_self, shape_tau, M, 
		 logRoDivOneMinRo, truncate_me, a_pareto_inv, m_minus2_div2, i, tau0);    
      // .. Build full B prior mat
      tau_i = tau_ij.row(i);
      priorMultiTau(priorMatB, smallPriorMat, tau_i, tau0, M, genes);
      // .. Calculate constants for gene i
      FullprecB   = priorMatB + eta(i)*C;            
      meanBNoPrec = eta(i)*(xt_plus1.row(i)-mu(i))*Y;
     
      // .. Update Gammas and Bs for row i 
      updateGammaAndB_row_i(B, gamma_ij, FullprecB, meanBNoPrec, logRosMinlogS, genes, M, i);
      
      // .. Correct identifiability mu, sum(YB)
      // *******************************************************************
      YB            = Y * trans(B.row(i));
      xhat_1(i)     = mean(YB);      
      xhat_0.col(i) = YB - xhat_1(i);    
      
      //.. Update mu(i)
      updateMu_Splines(mu, eta, eta_mu, B, mean_xt1, xhat_1, time_m, i);

      // .. Correct identifiability mu, sum(YB)
      // *******************************************************************
      mu(i) = mu(i) + xhat_1(i);
    }
    
    if(mcmc_iteration%10 == 0)
      estimateTime_AllowCancel(informTimeFlag, mcmc_iteration, samples, start_t);
    
    if((mcmc_iteration>burnIn) & (mcmc_iteration%thin == 0))
    {     
      // Timer? Estimated time?
      // Residuals?
      // Write variables to file
      writeMatToFile(TauFile, tau_ij);     
      writeToFileDouble(RhoFile, ro);      
      writeToFileVec(LambdaFile, eta);
      writeToFileVec(MuFile, mu);
      writeToFileInt(GammaFile, gamma_ij);   
      // .. Calculate value of functions
      calcIndivF(all_f, all_f_sqr, full_F_sqr, Y, B, genes, M, time_m);
      counter++;
    }    
  }   
  // .. Save mean and varaiance to file
  all_f   = all_f/counter;
  allf_file = ResultsFolder + "all_f"; 
  all_f.save(allf_file, raw_ascii);
  allf_sqr_file = ResultsFolder + "all_f_sqr"; 
  all_f_sqr = all_f_sqr/counter;
  all_f_sqr.save(allf_sqr_file, raw_ascii);
  full_F_sqr_file = ResultsFolder + "Full_F_sqr"; 
  full_F_sqr = full_F_sqr/counter;
  full_F_sqr.save(full_F_sqr_file, raw_ascii);  
  fclose(TauFile);
  fclose(RhoFile);
  fclose(LambdaFile);
  fclose(MuFile);
  fclose(GammaFile);
//   return 0;
}

