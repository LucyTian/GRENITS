// using namespace arma;
void paramFromVec_Splines(const arma::colvec& ParamVec_C, int& samples, int& burnIn, int& thin, double& c, 
		   double& d, double& a, double& b, double& truncate_me, double& tau0, int& M,
		   double& sigmaMu, double& a_pareto);
arma::mat despline(const arma::mat& data, const int& ndx, const int& deg);
arma::umat buildMapGammaBeta(int M, int genes);
void initMCMCvars_Splines(arma::colvec &mu, double &ro, arma::umat &gamma_ij, arma::mat &B, arma::colvec &eta, 
			  int genes, int conditions, arma::mat &tau_ij, int M);
void updateMu_Splines(arma::colvec& mu, const arma::colvec& eta, const double& eta_mu, const arma::mat& B, const arma::colvec& mean_xt1, 
		      const arma::rowvec& xhat_1, int time_m, int i);
void updateGammaAndB_row_i(arma::mat& B,   arma::umat& gamma_ij, const arma::mat& FullprecB, const arma::rowvec& meanBNoPrec, 
			     const arma::colvec & logRosMinlogS, int genes, int M, int i);

void openOutputFiles_Splines(string& ResultsFolder, FILE* &MuFile, FILE* &LambdaFile, FILE* &DesignFile, FILE* &GammaFile, 
		     FILE* &RhoFile, FILE* &TauFile, FILE* &All_fFile);
void fillBzerosUseGamma(arma::mat& B, const arma::umat& gamma_ij, int M);
arma::mat priprec(int M);
void updateTaus(arma::mat& tau_ij, arma::colvec& logRosMinlogS, const arma::mat& smallPriorMat, const arma::umat& gamma_ij, const arma::mat& B, 
		double shape_tau_self, double shape_tau, int M, double logRoDivOneMinRo, double truncate_me, double a_pareto_inv,
		double m_minus2_div2, int i, double tau0);
void priorMultiTau(arma::mat& priorWithTau, const arma::mat& smallPriorMat, const arma::rowvec& tau_i, double tau0, int M, int genes);
void calcIndivF(arma::mat& all_f, arma::mat& all_f_sqr, arma::mat& full_F_sqr, const arma::mat& Y, const arma::mat& B, int genes, int M, int time_m);
