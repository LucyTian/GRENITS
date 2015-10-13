void paramFromVec_Gauss(const colvec& ParamVec_C, int& samples, int& burnIn, int& thin, double& c, 
		        double& d, double& a, double& b, double& a_exp, double& b_exp, 
		        double& sigmaS, double& sigmaMu, int& numExp, int& thinDataSample);

void openOutputFiles_Gauss(string& ResultsFolder,    FILE* &Bfile,    FILE* &MuFile, 
			   FILE* &RhoFile, FILE* &LambFile, FILE* &GammaFile, FILE* &LambExpFile);

void readDataBioreps_ReturnAll(arma::cube& allData, arma::colvec& Xhat_sqr, arma::mat& Xhat, arma::mat& N_it, 
			       arma::mat& Y, int& genes, int& time_m, const arma::mat& dataMat, int numDataSets);
void reCalcYstats(arma::mat &C, arma::mat &Cplus , const arma::mat &Yt, const arma::mat &Ytplus1, const arma::mat &mu_mat);
void update_lambdaExp(arma::colvec &lambda_Exp, const arma::colvec &Xhat_sqr, 
		      const arma::mat &Xhat, const arma::mat &N_it, const arma::colvec &shape_lamExp, double b_exp, 
		      const arma::mat& Y);
void update_Y(arma::mat& Y, arma::mat& Yt, arma::mat& Ytplus1, const arma::mat Xhat, const arma::mat N_it, 
	      const arma::colvec& lambda_Exp, const arma::colvec& eta, int time_m, const arma::colvec& mu, 
	      const arma::mat& B, int genes);

// Student Functions
void openOutputFiles_Student(string& ResultsFolder,    FILE* &Bfile,    FILE* &MuFile, FILE* &RhoFile, 
						  FILE* &LambFile, FILE* &GammaFile, FILE* &LambExpFile, FILE* &DegFreedomFile);
void paramFromVec_Student(const colvec& ParamVec_C, int& samples, int& burnIn, int& thin, double& c, 
		   double& d, double& a, double& b, double& a_exp, double& b_exp, double& a_deg, 
		   double& b_deg,double& sigmaS, double& sigmaMu, int& numExp, int& thinDataSample);
		   
void initMCMCvars_Student(arma::colvec &mu, double &ro, arma::umat &gamma_ij, arma::mat &B, arma::colvec &eta, 
			  arma::colvec &lambda_Exp, arma::colvec &degFreedom, int genes, int conditions);
void update_weights_t(arma::cube &w_itr, arma::cube &YminX_sqr, const arma::mat &Y, const arma::cube &dataSet, 
		      const arma::colvec &degFreedom, const arma::mat &lambda_exp, int reps, int time_m);
void update_LambdaExp_t(arma::colvec &lambda_exp, const arma::cube &YminX_sqr,  const arma::cube &w_itr, 
			const arma::colvec &shape_lamExp, double b_exp);

void update_MH_DegFreedom_t(arma::colvec &degFreedom, arma::ucolvec &numAccepted, double a_deg, double b_deg, 
			    const arma::cube &w_itr, int genes, const arma::colvec &numSamples_i, 
			    double a_RW);
void update_Y_tDist(arma::mat& Y, arma::mat& Yt, arma::mat& Ytplus1, const arma::cube& dataSets, const arma::cube& w_itr, 
		    const arma::colvec& lambda_Exp, const arma::colvec& eta, int time_m, const arma::colvec& mu, 
		    const arma::mat& B, int genes);
