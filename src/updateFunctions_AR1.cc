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
// #include <google/profiler.h>


void MHStep(urowvec& gamma_Row,  double& logMVPDF_Old, const unsigned int& j, const mat& lambxCPlusS,   const rowvec& lambxCplusIdot, const double& sumLogs);
void calc_logMVPDF_withLinks(double& logMVPDF, const mat& lambxCPlusS , const rowvec& lambxCplusIdot, urowvec& gamma_Row);
void makeParametersRegression_i(   mat &lambxCPlusS, rowvec &lambxCplusIdot,   urowvec &links_on, 
				 urowvec &updateVec,   const umat &updateMe,    ucolvec &regsVec, 
				              int i,           const mat &C,    const mat &Cplus, 
				  const colvec &eta,   const umat &gamma_ij, double prior_prec_B,
			        ucolvec &numRegsVec);



void openOutputFiles_AR1(string& ResultsFolder, FILE* &Bfile, FILE* &MuFile, FILE* &RhoFile, FILE* &LambFile, FILE* &GammaFile)
{  
   // .. Declare strings
   string B_name, Mu_name, rho_name, lamb_name, gamma_name; 
   // .. Read numbers into it.
   B_name       = ResultsFolder + "B_mcmc"; 
   Mu_name      = ResultsFolder + "Mu_mcmc"; 
   rho_name     = ResultsFolder + "Rho_mcmc"; 
   lamb_name    = ResultsFolder + "Lambda_mcmc"; 
   gamma_name   = ResultsFolder + "Gamma_mcmc";   
   // .. Open files
   Bfile      = fopen(B_name.c_str(), "w");
   MuFile     = fopen(Mu_name.c_str(), "w");
   RhoFile    = fopen(rho_name.c_str(), "w");
   LambFile   = fopen(lamb_name.c_str(), "w");
   GammaFile  = fopen(gamma_name.c_str(), "w");
}

void paramFromVec_AR1(const colvec& ParamVec_C, int& samples, int& burnIn, int& thin, double& c, 
		   double& d, double& a, double& b, double& sigmaS, double& sigmaMu)
{  
   samples = ParamVec_C(0);
   burnIn = ParamVec_C(1);
   thin = ParamVec_C(2);
   c = ParamVec_C(3);
   d = ParamVec_C(4);
   sigmaS = ParamVec_C(5);
   a = ParamVec_C(6);
   b = ParamVec_C(7);
   sigmaMu = ParamVec_C(8);   
}
    


// Initialise MCMC variables
void initMCMCvars_AR1(colvec &mu, double &ro, umat &gamma_ij, mat &B, colvec &eta, int genes, int conditions)
{
  double r_min    = 0.0001;
  double r_max    = 0.2;
  double pb_min   = -1;
  double pb_max   = 1;
  double lamb_min = 0.1;
  double lamb_max = 1;

  // .. Set size of vars
  B.set_size(genes, genes);
  gamma_ij.set_size(genes, genes);
  eta.set_size(genes);
  mu.set_size(genes);

  // .. Init MCMC variables    
  ro =  Rf_runif(r_min, r_max);  
  RandomBernVec(gamma_ij.memptr(), ro, genes*genes);
  RandomUniformVec(mu.memptr(), pb_min, pb_max, genes);
  RandomUniformVec(B.memptr(), pb_min, pb_max, genes*genes);
  RandomUniformVec(eta.memptr(), lamb_min, lamb_max, genes);
}

// .. Sample from mu  -----------------------------------------------
void updateMu_AR1(colvec& mu, const colvec& eta, const double& eta_mu, const mat& B, const colvec& mean_xt1, const colvec& mean_xt,const unsigned int &  time_m)
{
    // Declare vars
    colvec aux_muVar, mu_sd, mu_mean;
    
    aux_muVar = 1/(1 + eta_mu / (eta*time_m));
    mu_sd     = sqrt(aux_muVar/(eta*time_m));
    mu_mean   = aux_muVar % (mean_xt1 - B*mean_xt);   
    
    // .. Sample
    for(unsigned int loop_var = 0; loop_var < mu.n_elem ;loop_var++)
    { 
      mu(loop_var) = Rf_rnorm(mu_mean(loop_var), mu_sd(loop_var));
    }
}

void updateCoeffAndGibbsVars(mat& B,   umat& gamma_ij, const colvec& eta, const mat& C, const mat& Cplus, const mat& precMatrix, 
			     const double & logRosMinlogS, const unsigned int & genes)
{
    // Declare vars
    mat            lambxCplusIdot, lambxCPlusSReduced, lambxCPlusS;
    rowvec        lambxCplusIdotReduced, lambxCplusIdotReduced_aux;
    urowvec                                               links_on;
    double                                            logMVPDF_Old;
    unsigned int                                                 i;
    int                                                          j;
    ucolvec                                          seq_rnd(genes);
    
    // .. Loop over rows i of matrix
    for(i = 0; i < genes;i++)
    {   
      // .. Calculate constants for gene i  
      lambxCPlusS    = eta(i)*C+precMatrix;
      lambxCplusIdot = eta(i)*Cplus.row(i);
      // .. Links in row i that are on
      links_on = gamma_ij.row(i); 
     
      // .. Calculate logMVPDF for this link config
      calc_logMVPDF_withLinks(logMVPDF_Old, lambxCPlusS, lambxCplusIdot, links_on);      
      
      // .. Random permutation of update index
      random_intSequence(seq_rnd);
      
      // .. Loop over elements of row i and sample from gamma
      for(int j_loop = 0; j_loop < genes; j_loop++)
      {
	// .. Get permuted index
	j = seq_rnd[j_loop];

	// .. Self interactions are not updated
	if (i != j)
	{
	  // .. MH update link i,j We reuse logMVPDF_Old (value 
	  // .. for current state of gamma_row)
	  MHStep(links_on, logMVPDF_Old, j, lambxCPlusS, lambxCplusIdot, logRosMinlogS);
	}
      }      
      
      // .. Update gamma matrix with new row
      gamma_ij.row(i) = links_on;  
      // .. Update coefficients
      updateCoefficients(B,i, links_on, lambxCPlusS, lambxCplusIdot);
   }

}


void updateCoeffAndGibbsVars_reg(mat& B,   umat& gamma_ij, const colvec& eta, const mat& C, const mat& Cplus, double prior_prec_B, 
			     const double & logRosMinlogS, const unsigned int & genes, umat& UpdateMe,
			     ucolvec &numRegsVec,      umat &regMat)
{
    // Declare vars
    mat                            lambxCPlusSReduced, lambxCPlusS;
    rowvec        lambxCplusIdotReduced, lambxCplusIdotReduced_aux;
    urowvec                                        links_on(genes);
    double                                            logMVPDF_Old;
    unsigned int                                                 i;
    int                                                          j;
    ucolvec                                                seq_rnd;
    rowvec                                          lambxCplusIdot;
    urowvec                                             updateVec ;
    ucolvec                                                regsVec;
    
    // .. Loop over rows i of matrix
    for(i = 0; i < genes;i++)
    {   
      // Get vector with indices of regulators
      getRegsVec(regsVec, numRegsVec, regMat, i);
      
      // .. Calculate parameters for iteration 
      makeParametersRegression_i( lambxCPlusS, lambxCplusIdot,   links_on, updateVec, UpdateMe,
				      regsVec,              i,          C,     Cplus,      eta,
				     gamma_ij,   prior_prec_B, numRegsVec);
      
      // .. Calculate logMVPDF for this link config
      calc_logMVPDF_withLinks(logMVPDF_Old, lambxCPlusS, lambxCplusIdot, links_on);      
      
      // .. Random permutation of update index
      seq_rnd.set_size(numRegsVec(i));
      random_intSequence(seq_rnd);
      
      // .. Loop over elements of row i and sample from gamma
      for(int j_loop = 0; j_loop < seq_rnd.n_elem; j_loop++)
      {
	// .. Get permuted index
	j = seq_rnd[j_loop];
	// .. Update non fixed links (self interaction is by default fixed)
	if (updateVec(j))
	{
	  // .. MH update link i,j We reuse logMVPDF_Old (value 
	  // .. for current state of gamma_row)
	  MHStep(links_on, logMVPDF_Old, j, lambxCPlusS, lambxCplusIdot, logRosMinlogS);
	}
      }       
      // .. Update gamma matrix with new row
      fillMatRowWithIndx_u(gamma_ij, links_on, i, regsVec);

      // .. Update coefficients
      updateCoefficients_reg(B, i, links_on, lambxCPlusS, lambxCplusIdot, regsVec);      
   }
}


// .. Metropolis-Hastings move to update a gamma_ij
// .. NOTE: Using log porbabilities. Always calculate move from
// .. link off to on. If move is inverese then multiply by -1
void MHStep(    urowvec& gamma_Row,  double& logMVPDF_Old,       const unsigned int& j, 
	    const mat& lambxCPlusS,   const rowvec& lambxCplusIdot_i, const double& sumLogs)
{
  // .. Vars
  unsigned int                                  gamma_old;  
  double    signOfK, alfa, log_aux, hastingsRatio; 
  double     logMVPDF_1, logMVPDF_0, logMVPDF_new;
  
  // .. Get current value of gamma_ij  
  gamma_old    = gamma_Row(j);
  // .. Change current value of gamma_ij (new proposed config)
  if (gamma_old){ gamma_Row(j) = 0; }else{ gamma_Row(j) = 1;}

  // .. Calculate logMVPDF for this link config
  calc_logMVPDF_withLinks(logMVPDF_new, lambxCPlusS, lambxCplusIdot_i, gamma_Row );

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

//   cout << logMVPDF_new << " " << logMVPDF_Old << endl;
  // .. Calculate hastings ratio for move from gamma_ij
  // .. and change sign if signOfK is -1
  hastingsRatio = sumLogs + 0.5 *(logMVPDF_1 - logMVPDF_0);
  hastingsRatio = signOfK * hastingsRatio;
  alfa          = min_two(0., hastingsRatio);
  log_aux       = log(unif_rand()); 
    
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
  } 
}


void calc_logMVPDF_withLinks(double& logMVPDF, const mat& lambxCPlusS , const rowvec& lambxCplusIdot_i, urowvec& gamma_Row)
{
  // .. Vars
  unsigned int                                            num_on;  
  mat                                         lambxCPlusSReduced;
  rowvec                                   lambxCplusIdotReduced;
  uvec                                                   regsVec;
  
  // .. Make sure there is at least one link
  regsVec = find(gamma_Row);
//   num_on   = accu(gamma_Row);
  if(regsVec.n_elem>0){	
      // .. Calculate logMVPDF for current state of gamma_row
//       subMatFromVector(lambxCPlusSReduced, lambxCPlusS, gamma_Row);
      subMatFromIndices(lambxCPlusSReduced, lambxCPlusS, regsVec);
//       subVectorFromVector(lambxCplusIdotReduced, lambxCplusIdot_i, gamma_Row);
      subVectorFromIndices(lambxCplusIdotReduced, lambxCplusIdot_i, regsVec);
      MHlogMVPDF(logMVPDF, lambxCPlusSReduced, lambxCplusIdotReduced);
  }
  else{
      logMVPDF = 0;
  }
}

void makeParametersRegression_i(   mat &lambxCPlusS, rowvec &lambxCplusIdot,   urowvec &links_on, 
				 urowvec &updateVec,   const umat &updateMe,    ucolvec &regsVec, 
				              int i,           const mat &C,    const mat &Cplus, 
				  const colvec &eta,   const umat &gamma_ij, double prior_prec_B,
			        ucolvec &numRegsVec)
{
//   ucolvec indx_regs;
  subMatFromIndices(lambxCPlusS, C, regsVec);

  // .. Calculate constants for gene i 
  // precision mat vals
  lambxCPlusS        = eta(i)*lambxCPlusS;
  lambxCPlusS.diag() = lambxCPlusS.diag()+ prior_prec_B; 
  
  subVectorFromIndx_MatRow(lambxCplusIdot, Cplus, i, regsVec);

  lambxCplusIdot     = eta(i)*lambxCplusIdot;

  subVectorFromIndx_MatRow_u(updateVec, updateMe, i, regsVec);
  subVectorFromIndx_MatRow_u(links_on, gamma_ij, i, regsVec);      
}

