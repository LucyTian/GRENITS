// using namespace arma;
void paramFromVec_AR1(const arma::colvec& ParamVec_C, int& samples, int& burnIn, int& thin, double& c, 
		   double& d, double& a, double& b, double& sigmaS, double& sigmaMu);

void initMCMCvars_AR1(arma::colvec &mu, double &ro, arma::umat &gamma_ij, arma::mat &B, arma::colvec &eta, int genes, int conditions);
void updateMu_AR1(arma::colvec& mu, const arma::colvec& eta, const double& eta_mu, const arma::mat& B, const arma::colvec& mean_xt1, const arma::colvec& mean_xt,const unsigned int &  time_m);
void updateCoeffAndGibbsVars(arma::mat& B,   arma::umat& gamma_ij, const arma::colvec& eta, const arma::mat& C, const arma::mat& Cplus, const arma::mat& precMatrix, 
			     const double & logRosMinlogS, const unsigned int & genes);

// void openOutputFiles(string& ResultsFolder, ofstream& Bfile, ofstream& pFile, ofstream& RhoFile, ofstream& LambFile, ofstream& Lamb_tauSqFile, ofstream& GammaFile);
void openOutputFiles_AR1(string& ResultsFolder, FILE* &Bfile, FILE* &MuFile, FILE* &RhoFile, FILE* &LambFile, FILE* &GammaFile);

