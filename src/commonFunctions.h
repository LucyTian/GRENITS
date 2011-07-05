// Not sure if this is necessary
#ifndef CTIME_HEADER
#include <ctime>
#define CTIME_HEADER
#endif
double min_two(double a, double b);
arma::colvec mvnrnd(arma::colvec &Mu, arma::mat& Sigma);
void RandomBernVec(unsigned int* array_rand, double ro, int total_elems);
void RandomUniformVec(double* array_rand, double min_val, double max_val, int total_elems);
void MHlogMVPDF(double& logMVPDF, const arma::mat& Sigma, const arma::rowvec& Mu);
void updateEta(arma::colvec& eta, const arma::mat& residuals, const double& shape_eta, const double& b);
arma::mat ScaleData(arma::mat & dataMat);
arma::mat loadAndScaleData(const char* filename);
void writeMatToFile(FILE* BFile, const arma::mat & B);
void writeToFileInt(FILE* GammaFile,const arma::umat& gamma_ij);
void writeToFileDouble(FILE* RhoFile, const double ro);
void writeToFileVec(FILE* vFile, const arma::colvec& v);
void random_intSequence(arma::ucolvec& seq);
void estimateTime_AllowCancel(arma::ucolvec& informTimeFlag_vec, int iteration_k, int samples, clock_t& start);
arma::umat is_NaN_mat(arma::mat &A);
void processFixedGammas(mat &Gamma_fixed, int &num_fixedON, int &free_gammas, umat &UpdateMe, umat &gamma_ij);

