#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
// #include <armadillo>
#include <RcppArmadillo.h>
using namespace arma;
using namespace std;
/*     
g++ -O3 -o NumParents MakeUncertaintyMat.cpp ../GRENITS/src/ReadMCMCFiles.cpp -larmadillo -lm -Wall
*/

void SetSizeVector(string &line, unsigned int& flagSetSize, vec &meanPost);
void MapMat2ReducedVector(mat &fixedMat, uvec &mapMatVec, vec &meanPost,  unsigned int &flagSetSize);
void FillNumParentsMat(mat &numParents, uvec &numParents_k);


// Based on:
// http://stackoverflow.com/questions/415515/how-can-i-read-and-manipulate-csv-file-data-in-c
void getPosteriorMeanFromFile(string &filename, vec &meanPost)
{
//    ifstream data("../Example_LinearNet_1.0/chain1/Gamma_mcmc");
//    ifstream data("../Runs/SYSMO_test/test1/chain1/Gamma_mcmc");
   ifstream data(filename.c_str());
   unsigned int          totalLines = 0;
   unsigned int         flagSetSize = 1;
   
//    vec                     meanPost;
   string                          line;
   vec::iterator      meanPost_iter;
      
   while(getline(data, line))
   {
      stringstream  lineStream(line);
      string        cell;
      // .. If size has not been set, set it 
      SetSizeVector(line, flagSetSize, meanPost);
      meanPost_iter = meanPost.begin();
      while(getline(lineStream, cell, ','))
      {
	  *meanPost_iter += atof(cell.c_str());
	  meanPost_iter++;
      }
      totalLines++;
   }
   meanPost = meanPost/totalLines;
//    meanPost.print_trans();
}


// Based on:
// http://stackoverflow.com/questions/415515/how-can-i-read-and-manipulate-csv-file-data-in-c
void getPosteriorMeanFromFile_withNumParents(string &filename, vec &meanPost, mat &numParents, int genes, mat &fixedMat)
{
   ifstream data(filename.c_str());
   unsigned int          totalLines = 0;
   unsigned int       flagSetSize_1 = 1;
   unsigned int       flagSetSize_2 = 1;
   uvec             numParents_k(genes);
   uvec 	              mapMatVec;
   string                          line;
   vec::iterator          meanPost_iter;
   uvec::iterator        mapMatVec_iter;
   double                       value_n;   
     
   numParents.zeros(genes, genes + 1);
   while(getline(data, line))
   {
      stringstream  lineStream(line);
      string                    cell;
      // .. If size has not been set, set it 
      SetSizeVector(line, flagSetSize_1, meanPost);
      MapMat2ReducedVector(fixedMat, mapMatVec, meanPost, flagSetSize_2);
      numParents_k.zeros();
      meanPost_iter  = meanPost.begin();
      mapMatVec_iter = mapMatVec.begin();      
      while(getline(lineStream, cell, ','))
      {
	  value_n         = atof(cell.c_str());
	  *meanPost_iter += value_n;
	  numParents_k[*mapMatVec_iter]+=value_n; 
	  meanPost_iter++;
	  mapMatVec_iter++;
      }
      totalLines++;
      FillNumParentsMat(numParents, numParents_k);
   }
//    cout << numParents_k<< endl<<endl;
//    cout << mapMatVec;
   meanPost   = meanPost/totalLines; 
   numParents = numParents/totalLines;
}


void SetSizeVector(string &line, unsigned int &flagSetSize, vec &meanPost)
{
  if(flagSetSize)
  {
    stringstream   lineStream(line);
    string                     cell;
    int                 numElem = 0;
    while(getline(lineStream, cell, ','))
      numElem++;
    meanPost.zeros(numElem);
    flagSetSize = 0;
  }
}

// Map position of flat reduced matrix to what regression element
// belongs to. Only do this once. Use meanPost to set size.
void MapMat2ReducedVector(mat &fixedMat, uvec &mapMatVec, vec &meanPost,  unsigned int &flagSetSize)
{
  if (flagSetSize)
  {
    mapMatVec.set_size(meanPost.n_elem);
    uvec::iterator mapMatVec_iter = mapMatVec.begin();
    unsigned int  col_j, row_i; 
    for(col_j = 0; col_j != fixedMat.n_cols; col_j++)
    {
      for(row_i = 0; row_i != fixedMat.n_cols; row_i++)
      {     
	if(isnan(fixedMat(row_i, col_j)))
	{
	  *mapMatVec_iter = row_i;
	  mapMatVec_iter++;	
	}
      }    
    }
    flagSetSize = 0;
  }
}


void FillNumParentsMat(mat &numParents, uvec &numParents_k)
{
  for(unsigned int row_i = 0; row_i != numParents.n_rows; row_i++)
      numParents(row_i, numParents_k[row_i]) ++ ;
}