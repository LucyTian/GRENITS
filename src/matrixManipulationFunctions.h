// using namespace arma;
// void subVector_FromIndices(arma::colvec& reducedVec,const arma::colvec& FullColVec, const arma::ucolvec& aux_blockMe);
void subVector_ShedIndices(arma::colvec& reducedVec,const arma::colvec& FullColVec, const arma::ucolvec& aux_blockMe);
void placeInVec(arma::rowvec& B_i_aux, const arma::colvec& B_block, const arma::ucolvec& aux_blockMe);
void getCAndC_offDiag(arma::mat& reduced_C, arma::mat& C_offDiag, const arma::mat& C, const arma::ucolvec& aux_blockMe);
void shedRow(arma::mat& Y , const arma::mat& X, const int& i);
void shedRowAndColumn(arma::mat& Y , const arma::mat& X, const int& i, const int& genes);
void subMatFromVector(arma::mat& reduced_C, const arma::mat& C, const arma::urowvec& logicVector);
void subMatFromIndices(arma::mat& reduced_C, const arma::mat& C, const arma::ucolvec& indxVals);
void subVectorFromVector(arma::rowvec& reduced_v, const arma::rowvec& v, const arma::urowvec& logicVector);
void subVectorFromIndices(arma::rowvec& reduced_v, const arma::rowvec& v, const arma::uvec& indxVector);
void subVectorFromVector_u(arma::urowvec& reduced_v, const arma::urowvec& v, const arma::urowvec& logicVector);
void fillMatRowWithVecAndZeros(arma::mat& B, const arma::colvec& b_new, const int&i, const arma::urowvec& logicVector);
void fillMatRowWithVecAndZeros_withIndex(arma::mat& B, const arma::colvec& b_new, const int&i, 
					 const arma::urowvec& logicVector, const arma::ucolvec& indxRegs);
arma::mat reorderMatColsFromVec(const arma::mat& B, const arma::colvec& newOrder);
arma::rowvec generate_seq(int start, int end);
arma::mat DiagnalBlockMat(const arma::mat& A, int num_repeats);
void DiagnalBlockMat2(arma::mat &BlockDiagMat, const arma::mat& A, int num_repeats);

void placeVecInVec_FromVec(arma::colvec& v, const arma::colvec& u, const arma::ucolvec& replaceMe);
void fillMatRowWithVec_u(arma::umat& B, const arma::urowvec& b_new, const int&i, const arma::urowvec& logicVector);
void fillMatRowWithIndx_u(arma::umat& B, const arma::urowvec& b_new, const int&i, const arma::ucolvec& indxVec);
void subVectorFromIndx_MatRow(arma::rowvec &outVec, const arma::mat& B, const int&i, const arma::ucolvec& indxRegs);
void subVectorFromIndx_MatRow_u(arma::urowvec &outVec, const arma::umat& B, const int&i, const arma::ucolvec& indxRegs);
void symmetriseMat(arma::mat &A);

void prod_Diag(double &result, arma::mat &A);
void modulus_ColVec(double &result, arma::vec &v);

// void subVectorFromIndx_MatRow_fromBases(arma::rowvec &outVec, const arma::mat& B, const int&i, const arma::ucolvec& indxRegs, int M)
void subMatFrom_RowIndices(arma::mat& reduced_C, const arma::mat& C, const arma::ucolvec& indxVals);
void subMatFrom_ColIndices(arma::mat& reduced_C, const arma::mat& C, const arma::ucolvec& indxVals);
