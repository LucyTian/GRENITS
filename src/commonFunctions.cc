#include <iostream>
#include <fstream>
#include <iomanip>
// #include <armadillo>
#include <RcppArmadillo.h>
#include <R_ext/Utils.h>
#include <cstdio>
#include "matrixManipulationFunctions.h"

using namespace std;
using namespace arma;
#include <string>
// Not sure if this is necessary
#ifndef CTIME_HEADER
#include <ctime>
#define CTIME_HEADER
#endif


// LAPACK function
extern "C" void arma_fortran(dtrtrs)(char *uplo, char *trans, char *diag, blas_int *n, blas_int *nrhs, double *a, 
				     blas_int *lda, double *b, blas_int * ldb, blas_int *info);  
void solve_Tri_LAPACK(mat& R, colvec& v);
void estimateRemainingTime(double& percent_done, double& time_left, int length, 
			   int iteration_k, clock_t& start);
umat is_NaN_mat(mat &A);

			   
// .. Sample from eta (inverse sigma square) ---------------------------
void updateEta(colvec& eta, const mat& residuals, const double& shape_eta, const double& b)
{
    // Declare vars
    colvec sum_squares, rate_eta, scale_eta;
    // .. Aux variables for eta  
    sum_squares = sum(square(residuals),1); 
    rate_eta    = b + 0.5 * sum_squares;
    scale_eta   = 1/rate_eta;
    
    // .. Sample
    for(unsigned int loop_var = 0; loop_var < rate_eta.n_elem ;loop_var++)
    {
      eta(loop_var) = Rf_rgamma(shape_eta, scale_eta(loop_var));
    }
}

// Read data, subtract mean of each row and divide by stdv
mat ScaleData(mat & dataMat)
{  
  // Declare Vars
  mat meanMat, stdMat, mean_vec, std_vec;
  // Mean and stdev
  mean_vec = mean(dataMat, 1);
  std_vec  = stddev(dataMat, 0, 1);
  // Make mat from vecs
  meanMat  = repmat(mean_vec, 1, dataMat.n_cols);  
  stdMat   = repmat(std_vec, 1, dataMat.n_cols);
  // Scale
  dataMat  = (dataMat-meanMat)/stdMat;
  return(dataMat);
}


// Read data, subtract mean of each row and divide by stdv
mat loadAndScaleData(const char* filename)
{  
  // Declare Vars
  mat dataMat, meanMat, stdMat, mean_vec, std_vec;
  // Load data matrix
  dataMat.load(filename);
  // Mean and stdev
  mean_vec = mean(dataMat, 1);
  std_vec  = stddev(dataMat, 0, 1);
  // Make mat from vecs
  meanMat  = repmat(mean_vec, 1, dataMat.n_cols);  
  stdMat   = repmat(std_vec, 1, dataMat.n_cols);
  // Scale
  dataMat  = (dataMat-meanMat)/stdMat;
  return(dataMat);
}


// // .. log of multivariate normal pdf as used in MH update
// void MHlogMVPDF(double& logMVPDF, const mat& Sigma, const rowvec& Mu)
// {
//   // .. Vars
//   mat          invR, R;
//   rowvec     vecDotMat;
//   
//   //.. Use choleski to speed up
//   R         = chol(Sigma);
// 
//   invR      = inv(R);
//   vecDotMat = Mu*invR;
//   logMVPDF  = 2*log(prod(invR.diag())) + dot(vecDotMat,vecDotMat);
// }

// // .. log of multivariate normal pdf as used in MH update
// void MHlogMVPDF(double& logMVPDF, const mat& Sigma, const rowvec& Mu)
// {
//   // .. Vars
//   mat                 invR, R;
//   colvec     vecDotMat, MuAux;
//   
//   //.. Use choleski to speed up
//   R         = chol(Sigma);
//   
//   MuAux = trans(Mu);
// //   solve(vecDotMat, trans(R), MuAux);
//   solve_Tri_LAPACK(R, MuAux);
//   logMVPDF  = -2*log(prod(R.diag())) + dot(MuAux,MuAux);
// }

// .. log of multivariate normal pdf as used in MH update
void MHlogMVPDF(double& logMVPDF, const mat& Sigma, const rowvec& Mu)
{
  // .. Vars
  mat                     invR, R;
  colvec         vecDotMat, MuAux;
  double   proDiagonal, modulusMu;
  //.. Use choleski to speed up
  R         = chol(Sigma);
  
  MuAux = trans(Mu);

  solve_Tri_LAPACK(R, MuAux);

  prod_Diag(proDiagonal, R);
  modulus_ColVec(modulusMu, MuAux);
  logMVPDF  = -2*log(proDiagonal) + modulusMu;

}

void solve_Tri_LAPACK(mat& R, colvec& v)
{
  blas_int           m = R.n_rows;
  char            upper_tri = 'U';
  char                   nu = 'N';
  char                  trs = 'T';
  blas_int               info_out;
  blas_int                 nrhs=1;

  arma_fortran(dtrtrs)( &upper_tri, &trs, &nu, &m, &nrhs, R.memptr(), 
                                     &m, v.memptr(), &m, &info_out);      
}

void RandomUniformVec(double* array_rand, double min_val, double max_val, int total_elems)
{
   // Generate random numbers
  for(int random_loop = 0; random_loop < total_elems; random_loop++)
  {
    array_rand[random_loop] = Rf_runif(min_val, max_val);
  }
}

void RandomBernVec(unsigned int* array_rand, double ro, int total_elems)
{
   // Generate random numbers
  for(int random_loop = 0; random_loop < total_elems; random_loop++)
  {
    array_rand[random_loop] = Rf_rbinom(1, ro);
  }  
}

void writeMatToFile(FILE* BFile, const mat & B)
{
  mat::const_iterator a = B.begin();
  mat::const_iterator b = B.end();
  b--;  
  for(mat::const_iterator i = a; i != b; ++i)
  {
      fprintf(BFile, "%4.3f,", *i);
  }
  fprintf(BFile, "%4.3f\n", *b);
}

void writeMatToFile_withIndx(FILE* BFile, const mat & B, uvec &UpdateIndx_Vec)
{
  uvec::iterator a = UpdateIndx_Vec.begin();
  uvec::iterator b = UpdateIndx_Vec.end();
  b--;  
  for(uvec::const_iterator i = a; i != b; ++i)
  {
      fprintf(BFile, "%4.3f,", B[*i]);
  }
  fprintf(BFile, "%4.3f\n", B[*b]);
}

void writeToFileDouble(FILE* RhoFile, const double ro)
{
    fprintf(RhoFile, "%4.3f\n", ro);
}

void writeToFileVec(FILE* vFile, const colvec& v)
{
  colvec::const_iterator a = v.begin();
  colvec::const_iterator b = v.end();
  b--;
  for(colvec::const_iterator i = a; i != b; ++i)
  {
    fprintf(vFile, "%4.3f,", *i);
  }
  fprintf(vFile, "%4.3f\n", *b);
}

void writeToFileInt(FILE* GammaFile, const umat& gamma_ij)
{
  umat::const_iterator a = gamma_ij.begin();
  umat::const_iterator b = gamma_ij.end();
  b--;  
  for(umat::const_iterator i = a; i != b; ++i)
  {
    fprintf(GammaFile, "%d,", *i);
  }
  fprintf(GammaFile, "%d\n", *b);
}

void writeToFileInt_withIndx(FILE* GammaFile, const umat& gamma_ij, uvec &UpdateIndx_Vec)
{
  uvec::iterator a = UpdateIndx_Vec.begin();
  uvec::iterator b = UpdateIndx_Vec.end();
  b--;  
  for(umat::const_iterator i = a; i != b; ++i)
  {
    fprintf(GammaFile, "%d,", gamma_ij[*i]);
  }
  fprintf(GammaFile, "%d\n", gamma_ij[*b]);
}


double min_two(double a, double b)
{
  if (a<b)
  {
    return(a);
  }else{
    return(b); 
  }
}

colvec mvnrnd(colvec &Mu, mat& Sigma)
{
  // Declare vars
  unsigned int            loop_var;
  mat                   chol_sigma;
  colvec         output(Mu.n_elem);
  rowvec unit_gaussians(Mu.n_elem);
  // Sample from unit gaussians  
  for(loop_var = 0; loop_var < Mu.n_elem;loop_var++)
  {    
    unit_gaussians[loop_var] = norm_rand();    
  }  
  
  chol_sigma = chol(Sigma);
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Z|ZZZZZZZZZZZZZZZZZZZZZZ
/*  if(chol_sigma.n_cols==0)
  {
    cout << "Size zero" <<endl;
    cout << "Num cols Sigma" << Sigma.n_cols <<endl;
//     Sigma.print("SigmaMat");
  }*/
  output = Mu + trans(unit_gaussians * chol_sigma) ;
  return(output);  
}

// .. From Durstenfeld (Algorithm 235: Random permutation)
void random_intSequence(ucolvec& seq)
{
  int r_int, old_val;
//   int a[seqSize];  
  for (int i = 0; i < seq.n_elem; i++)
  {
    seq[i] = i;
  }
  for (int i = seq.n_elem - 1; i!=0; --i)
  {
    r_int = floor(Rf_runif(0, i+1));
    old_val  = seq[r_int];
    seq[r_int] = seq[i]; 
    seq[i]     = old_val;
  }
}


void estimateTime_AllowCancel(ucolvec& informTimeFlag_vec, int iteration_k, int samples, 
			      clock_t& start)
{
  double percent_done, time_left;
  // .. Estimate run time 
  estimateRemainingTime(percent_done, time_left, samples, iteration_k, start);
  // .. Check to see if intial time estimate was given
  if(informTimeFlag_vec[0])
  {
    // .. If different from -1
    if(time_left > 0.)
    {
      Rcpp::Rcout << "Estimated runtime = " << time_left << " min" << endl;
      // Summary has already been given
      informTimeFlag_vec[0] = 0;
    }  
  }
  if(percent_done==25 or percent_done==50 or percent_done==75)
  {
    Rcpp::Rcout << percent_done << "% completed" << endl;
  }
  
  // .. Allow for user interrupt
  R_CheckUserInterrupt();  
}


void estimateRemainingTime(double& percent_done, double& time_left, int length, 
			        int iteration_k,    clock_t& start)
{
  double so_far, percent_speed, fraction_done;
  double so_far_sec;
//   long                   so_far_sec, time_now;
  clock_t                          so_far_clock, time_now; //so_far_sec, 
  
  // .. Percentage done 
  fraction_done = (double)iteration_k/(double)length;
  percent_done  = fraction_done*100;

  // .. Remaining time 
  time_now       = clock();
  so_far_clock   = time_now-start;
  so_far_sec     = (double)(time_now-start)/CLOCKS_PER_SEC;
  so_far         = so_far_sec/60.;
  percent_speed  = percent_done/so_far;
  
  if (so_far_sec > 1.)
  {
    time_left      = (100. - percent_done)/percent_speed;
  }
  else
  {
    time_left = -1.;
  }
}
// .. Take update information matrix. Elements are 0/1 for fixed off/on 
// .. and NaN for updatable link.  
void processFixedGammas(mat &Gamma_fixed,     int &num_fixedON,  int &free_gammas,  umat &UpdateMe, 
			  umat &gamma_ij,  ucolvec &numRegsVec,      umat &regIndx_Mat, 
			  uvec &UpdateIndx_Vec, uvec &flatRegsIndx_Vec)
{
  int indx_RegIndxMat;
  regIndx_Mat.zeros(Gamma_fixed.n_cols, Gamma_fixed.n_cols);
  umat   regMat;
  UpdateMe    = is_NaN_mat(Gamma_fixed);
  // Init regMat (position of regulators) as UpdateMe
  regMat      = UpdateMe;
  free_gammas = accu(UpdateMe);
//   Gamma_fixed.save("FixedGammaFile", raw_ascii); 
  num_fixedON = 0;
  // .. Place fixed values in gamma_ij and calculate how many are fixed on
  mat::iterator                         gamma_fixed_iter;
  mat::iterator        first_elem =  Gamma_fixed.begin();
  mat::iterator         last_elem =    Gamma_fixed.end();
  umat::iterator   iter_gamma_var =     gamma_ij.begin();
  umat::iterator      regMat_iter =       regMat.begin();
  
  for( gamma_fixed_iter = first_elem; gamma_fixed_iter != last_elem; 
        gamma_fixed_iter++, iter_gamma_var++, regMat_iter++)
  {
    if(!isnan(*gamma_fixed_iter))
    {
	*iter_gamma_var = *gamma_fixed_iter;
	num_fixedON += *gamma_fixed_iter;
	// If it is fixed on it's a regulator
	if(*gamma_fixed_iter==1)
	  *regMat_iter = 1;
    }
  } 
  // .. Number of regulators for each regression
  numRegsVec = sum(regMat,1);
  
  // .. Make Reg indx mat (fill with zeros place indices in col wise order)
  for(int i = 0; i!=Gamma_fixed.n_cols;i++)
  {
    indx_RegIndxMat=0;
    for(int j = 0; j!=Gamma_fixed.n_cols;j++)
    {
      if(regMat(i,j))
      {  
	regIndx_Mat(indx_RegIndxMat,i) = j;
	indx_RegIndxMat++;
      }
    }
  }
  UpdateIndx_Vec   = find(UpdateMe);
  flatRegsIndx_Vec = find(regMat); 
}

umat is_NaN_mat(mat &A)
{
  umat isNAN_vals = zeros<umat>(A.n_rows, A.n_cols);
  
  mat::iterator                              A_iter;
  mat::iterator     first_elem =          A.begin();
  mat::iterator      last_elem =            A.end();
  umat::iterator      iter_NAN = isNAN_vals.begin();
  
  for( A_iter = first_elem; A_iter != last_elem; A_iter++, iter_NAN++)
  {
    if(isnan(*A_iter))
    {
	*iter_NAN = 1;
    }
  } 
  return(isNAN_vals);
}

void getRegsVec(ucolvec &regsVec, ucolvec &numRegs, umat &regMat, unsigned int i)
{  
   unsigned int start_pos;
   // Calc flat mat index for col i row 0
   start_pos = i*regMat.n_rows;
   regsVec.set_size(numRegs(i));
   for(unsigned int loop_j = 0; loop_j< numRegs(i); loop_j++)
    regsVec[loop_j] = regMat[start_pos + loop_j];  
}


void getRegsVecBases(ucolvec &regsVec, ucolvec &numRegs, umat &regMat, unsigned int i, int M)
{  
   unsigned int start_pos;
   // Calc flat mat index for col i row 0
   start_pos = i*regMat.n_rows;
   regsVec.set_size(numRegs(i)*M);
   for(unsigned int loop_j = 0; loop_j< M*numRegs(i); loop_j++)
    regsVec[loop_j] = regMat[start_pos + loop_j];  
}

// .. C and X are NOT passed by reference
void updateCoefficients_reg(   mat& B, const int& i, const urowvec& links_on, const mat& lambxCPlusS, 
			       const rowvec& lambxCplusIdot_i, const ucolvec& indxRegs)
{
  // .. Declare vars
  mat                     covarianceMatrixB,  lambxCPlusSReduced;
  rowvec        lambxCplusIdotReduced, lambxCplusIdotReduced_aux;
  unsigned int                                            num_on;
  colvec                                         aux_mu_B, b_new;
  uvec                                                   regsVec;
  
  // .. Make sure there is at least one link
  regsVec = find(links_on);
  if(regsVec.n_elem>0){	
      // .. Calculate logMVPDF for current state of gamma_row
      subMatFromIndices(lambxCPlusSReduced, lambxCPlusS, regsVec);
      subVectorFromIndices(lambxCplusIdotReduced, lambxCplusIdot_i, regsVec);

//   // .. Before updating links, check if any need updating
//   num_on   = accu(links_on);
//   if (num_on>0)
//   {
//       // .. Aux variables for B    
//       subMatFromVector(lambxCPlusSReduced, lambxCPlusS, links_on);
//       subVectorFromVector(lambxCplusIdotReduced, lambxCplusIdot_i, links_on);


      // .. B Covariance matrix
      covarianceMatrixB = inv(lambxCPlusSReduced);
      // .. Correct precision => symmetry problem
//       covarianceMatrixB = (covarianceMatrixB + trans(covarianceMatrixB))/2;
      symmetriseMat(covarianceMatrixB);  
      // .. Update B for links that are ON
      aux_mu_B  = covarianceMatrixB*trans(lambxCplusIdotReduced);
      b_new     = mvnrnd(aux_mu_B, covarianceMatrixB);
  }
  // .. Update B row with new values and zeros
  fillMatRowWithVecAndZeros_withIndex(B, b_new, i, links_on, indxRegs);
}



// .. C and X are NOT passed by reference
void updateCoefficients(                mat& B,              const int& i, const urowvec& links_on, 
			const mat& lambxCPlusS, const rowvec& lambxCplusIdot_i)
{
  // .. Declare vars
  mat                     covarianceMatrixB,  lambxCPlusSReduced;
  rowvec        lambxCplusIdotReduced, lambxCplusIdotReduced_aux;
  unsigned int                                            num_on;
  colvec                                         aux_mu_B, b_new;
  uvec                                                   regsVec;
  
  
  // .. Make sure there is at least one link
  regsVec = find(links_on);
  if(regsVec.n_elem>0){	
      // .. Calculate logMVPDF for current state of gamma_row
      subMatFromIndices(lambxCPlusSReduced, lambxCPlusS, regsVec);
      subVectorFromIndices(lambxCplusIdotReduced, lambxCplusIdot_i, regsVec);

//   // .. Before updating links, check if any need updating
//   num_on   = accu(links_on);
//   if (num_on>0)
//   {
//       // .. Aux variables for B    
//       subMatFromVector(lambxCPlusSReduced, lambxCPlusS, links_on);
//       subVectorFromVector(lambxCplusIdotReduced, lambxCplusIdot_i, links_on);

      // .. B Covariance matrix
      covarianceMatrixB = inv(lambxCPlusSReduced);
      // .. Correct precision => symmetry problem
//       covarianceMatrixB = (covarianceMatrixB + trans(covarianceMatrixB))/2;
      symmetriseMat(covarianceMatrixB);  

      // .. Update B for links that are ON
      aux_mu_B  = covarianceMatrixB*trans(lambxCplusIdotReduced);
      b_new     = mvnrnd(aux_mu_B, covarianceMatrixB);
  }
  // .. Update B row with new values and zeros
  fillMatRowWithVecAndZeros(B, b_new, i, links_on);
}


