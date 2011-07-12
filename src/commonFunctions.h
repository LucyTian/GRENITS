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
void writeMatToFile_withIndx(FILE* BFile, const arma::mat & B, arma::uvec &UpdateIndx_Vec);
void writeToFileInt(FILE* GammaFile,const arma::umat& gamma_ij);
void writeToFileInt_withIndx(FILE* GammaFile,const arma::umat& gamma_ij, arma::uvec &UpdateIndx_Vec);
void writeToFileDouble(FILE* RhoFile, const double ro);
void writeToFileVec(FILE* vFile, const arma::colvec& v);
void random_intSequence(arma::ucolvec& seq);
void estimateTime_AllowCancel(arma::ucolvec& informTimeFlag_vec, int iteration_k, int samples, clock_t& start);
arma::umat is_NaN_mat(arma::mat &A);
void processFixedGammas( arma::mat &Gamma_fixed,     int &num_fixedON,  int &free_gammas,  arma::umat &UpdateMe, 
			arma::umat    &gamma_ij, arma::ucolvec  &numRegsVec, arma::umat &regMat,
		        arma::uvec &UpdateIndx_Vec, arma::uvec &flatRegsIndx_Vec);
void getRegsVec(arma::ucolvec &regsVec, arma::ucolvec &numRegs, arma::umat &regMat, unsigned int i);
void getRegsVecBases(arma::ucolvec &regsVec, arma::ucolvec &numRegs, arma::umat &regMat, unsigned int i, int M);

void updateCoefficients(                arma::mat& B,                       const int& i, const arma::urowvec& links_on, 
			const arma::mat& lambxCPlusS, const arma::rowvec& lambxCplusIdot);

void updateCoefficients_reg(                        arma::mat& B,                         const int& i, const arma::urowvec& links_on, 
				    const arma::mat& lambxCPlusS, const arma::rowvec& lambxCplusIdot_i, const arma::ucolvec& indxRegs);

