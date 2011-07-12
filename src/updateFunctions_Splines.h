// using namespace arma;
void paramFromVec_Splines(const arma::colvec& ParamVec_C, int& samples, int& burnIn, int& thin, double& c, 
		   double& d, double& a, double& b, double& truncate_me, double& tau0, int& M,
		   double& sigmaMu, double& a_pareto);
arma::mat despline(const arma::mat& data, const int& ndx, const int& deg);
arma::umat buildMapGammaBeta(int M, int genes);
void initMCMCvars_Splines(arma::colvec &mu, double &ro, arma::umat &gamma_ij, arma::mat &B, arma::colvec &eta, 
			  int genes, int conditions, arma::mat &tau_ij, int M);
void updateMu_Splines(arma::colvec& mu, const arma::colvec& eta, const double& eta_mu, const arma::colvec& mean_xt1, 
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
void updateTaus_reg(arma::mat& tau_ij, arma::colvec& logRosMinlogS, const arma::mat& smallPriorMat, arma::urowvec& links_on, const arma::mat& B, 
		double shape_tau_self, double shape_tau, int M, double logRoDivOneMinRo, double truncate_me, double a_pareto_inv,
		double m_minus2_div2, int i, double tau0, arma::uvec regsVec);
void priorMultiTau(arma::mat& priorWithTau, const arma::mat& smallPriorMat, const arma::rowvec& tau_i, double tau0, int M, int genes);
void calcIndivF(arma::mat& all_f, arma::mat& all_f_sqr, arma::mat& full_F_sqr, const arma::mat& Y, const arma::mat& B, int genes, int M, int time_m);
void fixedBasesFromFixedRegs(arma::umat &fixedBases, arma::umat &regMat, arma::ucolvec &numRegs, int M);
void makeParametersSplinesRegression_i(   arma::mat &FullprecB, arma::rowvec &meanBNoPrec,  arma::urowvec &updateVec,   
			            const arma::umat &updateMe,    arma::ucolvec &regsVec,               int i,           
				            const arma::mat &C,          arma::mat &Y_red,   const arma::colvec &eta,   
			              arma::mat &smallPriorMat,         arma::mat& tau_ij, arma::ucolvec &numRegsVec, 
				   arma::ucolvec &regsBasesVec,             arma::vec &mu,               double tau0,                     
				                         int M,              arma::mat &Y,      arma::mat  &xt_plus1);
void updateGammaAndB_row_i_reg(arma::mat& B,   arma::umat& gamma_ij, const arma::mat& FullprecB, const arma::rowvec& meanBNoPrec, 
			     const arma::colvec & logRosMinlogS, int genes, int M, int i, arma::urowvec &links_on,
			     const arma::ucolvec &regsBasesVec,  const arma::urowvec &updateVec, const arma::ucolvec &numRegsVec,
			     const arma::uvec &regsVec);
// void updateGammaAndB_row_i_reg(mat& B,   umat& gamma_ij, const mat& FullprecB, const rowvec& meanBNoPrec, 
// 			     const colvec & logRosMinlogS, int genes, int M, int i, urowvec &links_on,
// 			     const ucolvec &regsBasesVec,  const urowvec &updateVec, const ucolvec &numRegsVec,
// 			     const uvec &regsVec)
