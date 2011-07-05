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

// .. log of multivariate normal pdf as used in MH update
void MHlogMVPDF(double& logMVPDF, const mat& Sigma, const rowvec& Mu)
{
  // .. Vars
  mat                 invR, R;
  colvec     vecDotMat, MuAux;
  
  //.. Use choleski to speed up
  R         = chol(Sigma);
  
  MuAux = trans(Mu);
//   solve(vecDotMat, trans(R), MuAux);
  solve_Tri_LAPACK(R, MuAux);
  logMVPDF  = -2*log(prod(R.diag())) + dot(MuAux,MuAux);
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
    unit_gaussians(loop_var) = norm_rand();    
  }  
  
  chol_sigma = chol(Sigma);
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
      cout << "Estimated runtime = " << time_left << " min" << endl;
      // Summary has already been given
      informTimeFlag_vec[0] = 0;
    }  
  }
  if(percent_done==25 or percent_done==50 or percent_done==75)
  {
    cout << percent_done << "% completed" << endl;
  }
  
  // .. Allow for user interrupt
  R_CheckUserInterrupt();  
}


void estimateRemainingTime(double& percent_done, double& time_left, int length, 
			   int iteration_k, clock_t& start)
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

void processFixedGammas(mat &Gamma_fixed, int &num_fixedON, int &free_gammas, umat &UpdateMe, umat &gamma_ij)
{
  UpdateMe    = is_NaN_mat(Gamma_fixed);
  free_gammas = accu(UpdateMe);
  Gamma_fixed.save("FixedGammaFile", raw_ascii); 
  num_fixedON = 0;
  // .. Place fixed values in gamma_ij and calculate how many are fixed on
  mat::iterator                         gamma_fixed_iter;
  mat::iterator        first_elem =  Gamma_fixed.begin();
  mat::iterator         last_elem =    Gamma_fixed.end();
  umat::iterator   iter_gamma_var =     gamma_ij.begin();
  
  for( gamma_fixed_iter = first_elem; gamma_fixed_iter != last_elem; gamma_fixed_iter++, iter_gamma_var++)
  {
    if(!isnan(*gamma_fixed_iter))
    {
	*iter_gamma_var = *gamma_fixed_iter;
	num_fixedON += *gamma_fixed_iter;
    }
  } 
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


