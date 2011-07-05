#include <iostream>
// #include <armadillo>
#include <RcppArmadillo.h>

using namespace  std;
using namespace arma;

// .. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// .. 
// ..    Take matrix, vector, row index and [0,1] vector. Return Matrix with row 
// ..    replaced by vector using indeces given by [0,1] vector, zeros elswhere in
// ..    that row
// .. 
// .. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void fillMatRowWithVecAndZeros(mat& B, const colvec& b_new, const int&i, const urowvec& logicVector)
{
  // .. Iterator for new vals
  colvec::const_iterator b_new_iter = b_new.begin();
  unsigned int row_elem, col_lengh, indx_row_elem;
  col_lengh = B.n_rows;
  for( row_elem = 0; row_elem <logicVector.n_elem; row_elem++)
  {
    // .. Calculate index of row element
    indx_row_elem = i + row_elem*col_lengh;
    // .. If there is a non zero value insert next value
    if (logicVector[row_elem])
    {
      B[indx_row_elem]   = *b_new_iter;
      b_new_iter++;
    // .. Otherwise insert a zero
    }else{
      B[indx_row_elem] = 0;
    }
  }  
}


// .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// .. 
// ..    Take Vector and [0,1] vector. Return subvector defined by [0,1] vector
// .. 
// .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void subVectorFromVector(rowvec& reduced_v, const rowvec& v, const urowvec& logicVector)
{
  // .. Vars
  ucolvec indxVals, col_Indx;
  
  // .. Get indexes of wanted elements
  indxVals = find(logicVector);

  // .. Set dimensions of submatrix
  reduced_v.set_size(indxVals.n_elem);
  
  // .. Iterators for vector
  ucolvec::iterator    first_row_iter = indxVals.begin();  
  ucolvec::iterator     last_row_iter =   indxVals.end();
  ucolvec::iterator                             row_iter;
  
  // .. Iterator for vector to fill
  rowvec::iterator       reduced_v_iter = reduced_v.begin();
  
  for( row_iter = first_row_iter; row_iter != last_row_iter; row_iter++)
  {
    // .. Use iterators to calc (flat) index of C
    *reduced_v_iter = v[*row_iter];
    reduced_v_iter++;
  } 
}

// .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// .. 
// ..    Take matrix and [0,1] vector. Return submatrix defined by vector
// .. 
// .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void subMatFromVector(mat& reduced_C, const mat& C, const urowvec& logicVector)
{
  ucolvec           indxVals,  col_Indx;
  unsigned int                 genes = C.n_cols;
  
  // .. Get indexes of wanted elements
  indxVals = find(logicVector);
  // .. Get first row elements of each wanted col
  col_Indx =    genes*indxVals;

  // .. Set dimensions of submatrix
  reduced_C.set_size(indxVals.n_elem, indxVals.n_elem);
  
  // .. Iterators for first elemnts of each col (full mat)
  ucolvec::iterator    first_elem_iter = col_Indx.begin();  
  ucolvec::iterator     last_elem_iter =   col_Indx.end();
  ucolvec::iterator                             cols_iter;

  // .. Iterators for row vals (vector)
  ucolvec::iterator    first_row_iter = indxVals.begin();  
  ucolvec::iterator     last_row_iter =   indxVals.end();
  ucolvec::iterator                             row_iter;
  
  // .. Iterator for matrix to fill
  mat::iterator       reduced_C_iter = reduced_C.begin();
  
  for( cols_iter = first_elem_iter; cols_iter != last_elem_iter; cols_iter++)
  {
    for( row_iter = first_row_iter; row_iter != last_row_iter; row_iter++)
    {
      // .. Use iterators to calc (flat) index of C
      *reduced_C_iter = C[*cols_iter + *row_iter];
      reduced_C_iter++;
    } 
  }  
}


// .. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// .. 
// ..  Take output-array of correct size, donor array and indices colvec.
// ..  Return donor array minus elements given by indices 
// .. 
// .. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void subVector_ShedIndices(colvec& reducedVec,const colvec& FullColVec, const ucolvec& aux_blockMe)
{
  // .. Iterators for full
  colvec::const_iterator         Begin_Iter = FullColVec.begin();  
  colvec::const_iterator           End_Iter = FullColVec.end();
  // .. Output matrix iterator
  colvec::iterator                     y_it = reducedVec.begin();
  // .. Index to skip
  ucolvec::const_iterator  skipMeIndex_Iter = aux_blockMe.begin();
  ucolvec::const_iterator   lastSkipMe_Iter = aux_blockMe.end() - 1;
  
  // .. Loop to fill vector, skip unwanted values using skipMeIndex iterator
  for(mat::const_iterator elemN = Begin_Iter; elemN != End_Iter; ++elemN )
  {
      if(*elemN != FullColVec(*skipMeIndex_Iter))
      {      
	*y_it = *elemN;
	y_it++;
      }
      else
      {
	if (skipMeIndex_Iter!=lastSkipMe_Iter){skipMeIndex_Iter++;}
      }
  }
}

// .. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// .. 
// ..  Return array with elements from donor array given by indices 
// .. 
// .. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void subVector_FromIndices(colvec& reducedVec,const colvec& FullColVec, const ucolvec& aux_blockMe)
{
  // .. Iterators for full
  ucolvec::const_iterator          Begin_Iter = aux_blockMe.begin();  
  ucolvec::const_iterator            End_Iter = aux_blockMe.end();
  // .. Output matrix iterator
  colvec::iterator                       y_it = reducedVec.begin();
   // .. Full matrix iterator
  colvec::const_iterator     iterZero_FullVec = FullColVec.begin();
 
   // .. Loop to fill vector
  for(ucolvec::const_iterator AddMeIndex_Iter = Begin_Iter; AddMeIndex_Iter != End_Iter; ++y_it, ++AddMeIndex_Iter)
  {
	*y_it = *(iterZero_FullVec + *AddMeIndex_Iter);
  }
}

// .. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// .. 
// ..  Return array with elements replaced from donor array given by indices 
// .. 
// .. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void placeInVec(rowvec& B_i_aux, const colvec& B_block, const ucolvec& aux_blockMe)
{
  // .. Iterators for block updated index (with res[pect to full vector) and value
  ucolvec::const_iterator  iter_col_i       = aux_blockMe.begin();  
  ucolvec::const_iterator  last_elem_i      = aux_blockMe.end();  
  colvec::const_iterator   iter_B_block     = B_block.begin();
  rowvec::iterator         iterZero_FullVec = B_i_aux.begin();  
  
  // .. Place updated values in full vector
  for(ucolvec::const_iterator iter_indx = iter_col_i; iter_indx != last_elem_i; ++iter_indx,  ++iter_B_block)
  {
    *(iterZero_FullVec + *iter_indx) = *iter_B_block;
  } 
}


// .. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// .. 
// ..  Take matrix and indices. Return two matrices:
// ..  1 matrix with all rows using indices and col using indices 
// ..  2 matrix with all rows using indices and col excluding indices 
// .. 
// .. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void getCAndC_offDiag(mat& reduced_C, mat& C_offDiag, const mat& C, const ucolvec& aux_blockMe)
{
  unsigned int                            indxLastBlockMe =  aux_blockMe.n_elem - 1;
  ucolvec     aux_jumpVec(indxLastBlockMe + 1), jumpColvec(indxLastBlockMe);
  unsigned int                                                              nextCol;
  unsigned int                                                     genes = C.n_cols;
  int                                                        selfBlock_Flag;
//   urowvec                                                    unit_elem(1);
  
  // .. Build jump iterator to go from wanted index to wanted index
  jumpColvec      = ucolvec(aux_blockMe.rows(1, indxLastBlockMe)) - ucolvec(aux_blockMe.rows(0, indxLastBlockMe - 1));
  aux_jumpVec(0)  = aux_blockMe(0);
  aux_jumpVec.rows(1, indxLastBlockMe) = jumpColvec;
    
  // .. Iterators for block updated index (with res[pect to full vector) and value
  ucolvec::iterator  first_elem_iter = aux_jumpVec.begin();  
  ucolvec::iterator   last_elem_iter = aux_jumpVec.end();
  ucolvec::iterator    jumpElem_iter; 
  
  mat::const_iterator       posFullMat_iter = C.begin();

  mat::iterator              reduced_C_iter = reduced_C.begin();
  mat::iterator              C_offDiag_iter = C_offDiag.begin();
  double                      selfBlock_col = C(0, aux_blockMe(0)); 

  // .. Number of elements to advance to jump from last row in col to first rwo in next col
  nextCol = genes - aux_blockMe(indxLastBlockMe); 
  
  unsigned int colJump_elem_idx = 0;
  // .. Past last element
  unsigned int last_elem_jump = jumpColvec.n_elem;
  for(unsigned int loop_cols = 0; loop_cols < genes; loop_cols++)
  {
    // .. Check if we are in self block col -----------------
    if(*posFullMat_iter == selfBlock_col)
    { 
      // .. Set flag to True 
      selfBlock_Flag = 1;
      // .. Move selfBlock_elem to next (self block) col
//       cout << colJump_elem_idx;      
      if(colJump_elem_idx!=last_elem_jump){
	selfBlock_col = *(posFullMat_iter + jumpColvec(colJump_elem_idx)*genes);
	colJump_elem_idx++;
      }
    }else{
      // .. Set flag to false
      selfBlock_Flag = 0;
    }
    // -----------------------------------------------------
    // .. Loop over row elements ===================================================
    for(jumpElem_iter = first_elem_iter; jumpElem_iter != last_elem_iter; ++jumpElem_iter)
    {
      // .. Move to new pos
      posFullMat_iter += *jumpElem_iter; 
      if(selfBlock_Flag)
      {
	// .. Place value reduced C matrix 
	*reduced_C_iter = *posFullMat_iter;    
	reduced_C_iter++;
      }else{
	// .. Place value in off-diag matrix 
	*C_offDiag_iter = *posFullMat_iter;    
	C_offDiag_iter++;
      }    
    } 
    // .. ==========================================================================
    // .. Place iterator at first wanted row of next col
    posFullMat_iter += nextCol;
  }
  
}


// .. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// .. 
// ..  Return mat with row i removed
// .. 
// .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void shedRow(mat &Y , const mat &X, const int& i)
{  
  mat::const_iterator      a = X.begin();
  mat::const_iterator        b = X.end();
  mat::iterator         y_it = Y.begin();
  rowvec                row_i = X.row(i);
  unsigned int                      notThisOne=0;
  unsigned int       last_elems = row_i.n_elem-1;
  for(mat::const_iterator elemN = a; elemN != b; ++elemN)
  {
//       cout << notThisOne << ":";
      if(*elemN!=row_i(notThisOne))
      {      
	*y_it = *elemN;
	y_it++;
      }
      else
      {
	if (notThisOne != last_elems){notThisOne++;}
      }	
  }
}

// .. +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// .. 
// ..  Return mat with row i and col i removed
// .. 
// .. ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void shedRowAndColumn(mat &Y , const mat &X, const int& i, const int& genes)
{ 
  // .. General vars
  rowvec                row_i =       X.row(i);
  unsigned int              notThisOne =              0;
  unsigned int             last_elems = row_i.n_elem-1;  
  // .. First block
  mat::const_iterator      a_firstBlock = X.begin();  
  mat::const_iterator      b_firstBlock = a_firstBlock + genes * i;      
  // .. Second block
  mat::const_iterator      a_secondBlock = b_firstBlock + genes;    
  mat::const_iterator      b_secondBlock = X.end();
  // .. Output matrix
  mat::iterator            y_it = Y.begin();
  
  for(mat::const_iterator elemN = a_firstBlock; elemN != b_firstBlock; ++elemN)
  {
      if(*elemN!=row_i(notThisOne))
      {      
	*y_it = *elemN;
	y_it++;
      }
      else
      {
	notThisOne++;
      }	
  }
  // .. Skip column i (and shift row index one forward)
  notThisOne++;
  
  for(mat::const_iterator elemN = a_secondBlock; elemN != b_secondBlock; ++elemN)
  {
      if(*elemN!=row_i(notThisOne))
      {      
	*y_it = *elemN;
	y_it++;
      }
      else
      {
	if (notThisOne!=last_elems){notThisOne++;}
      }	
  }
}

mat reorderMatColsFromVec(const mat& B, const colvec& newOrder)
{
  // .. Mat to return
  mat return_mat(B.n_rows, B.n_cols);
  // .. Order vector iterator
  colvec::const_iterator order_iter = newOrder.begin();
  // .. Fill matrix
  for(unsigned int k = 0; k < B.n_cols; k++)
  {
    return_mat.col(k) = B.col(*order_iter);
    order_iter++;
  }
  return(return_mat);
}

rowvec generate_seq(int start, int end)
{
  rowvec              out_seq;
  int     seq_elem, final_dim;
  
  final_dim = end - start + 1;
  out_seq.set_size(final_dim);
  seq_elem = start;
  for(int loop_var = 0; loop_var < final_dim; loop_var++)
  {
    out_seq[loop_var] = (double) seq_elem;
    seq_elem++;
  }
  return(out_seq);
}

// .. Create block diagonal matrix with A repeated num_repeats times
mat DiagnalBlockMat(const mat& A, int num_repeats)
{
  // .. Vars
  double                                        num_elems_A = A.n_cols;
  mat         BlockDiagMat(A.n_cols*num_repeats, A.n_cols*num_repeats);  
  int                                             start_indx, end_indx;
  
  // .. Fill mat
  for(int block_i = 0; block_i < num_repeats;  block_i++) 
  {
    start_indx = block_i*num_elems_A;
    end_indx    = (block_i + 1)*num_elems_A - 1;
    BlockDiagMat.submat(start_indx, start_indx, end_indx, end_indx) = A; 
  }    
  return(BlockDiagMat);
}


void placeVecInVec_FromVec(arma::colvec& to_vec, const arma::colvec& from_vec, const arma::ucolvec& replaceMe)
{
  ucolvec replaceMe_indx = find(replaceMe);
  // .. Iterators for block updated index (with res[pect to full vector) and value
  ucolvec::const_iterator  start_elem_i = replaceMe_indx.begin();  
  ucolvec::const_iterator  last_elem_i  = replaceMe_indx.end();  
  // .. Place updated values in full vector
  for(ucolvec::const_iterator iter_indx = start_elem_i; iter_indx != last_elem_i; ++iter_indx)
      to_vec[*iter_indx] = from_vec[*iter_indx];
}


// void fillWithOnesAndZeros(ucolvec& v, map, const urowvec& logicVector)
// {
//   // .. Iterators for block updated index (with res[pect to full vector) and value
//   urowvec::const_iterator  iter_elem_i      = logicVector.begin();  
//   urowvec::const_iterator  last_elem_i      = logicVector.end(); 
//   
//   ucolvec::const_iterator  v_iter;   
//   
//   // .. Place updated values in full vector
//   for(logicVec_iter = iter_elem_i; logicVec_iter != last_elem_i; ++logicVec_iter, ++v_iter)
//   {
//     if(*logicVec_iter){
//       *v_iter = 1;
//     }else{
//       *v_iter = 0;
//     } 
// }

