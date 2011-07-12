// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// fastLm.cpp: Rcpp/Armadillo glue example of a simple lm() alternative
//
// Copyright (C)  2010 Dirk Eddelbuettel, Romain Francois and Douglas Bates
//
// This file is part of RcppArmadillo.
//
// RcppArmadillo is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppArmadillo is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppArmadillo.h>
#include "AR1_biocond.h"
#include "PSplines_biocond.h"
#include "AR1_Gauss_biocond.h"
#include "AR1_Student_biocond.h"
#include "ReadMCMCFiles.h"
#include <iostream>
#include <cstdio>
#include <string>


// #include <google/profiler.h>


extern "C" SEXP callAR1(SEXP Xs, SEXP ResFolder, SEXP ParamVec, SEXP FixedLinks) {
  try {
	// .. Transform parameters to c format ===================================================	
	// .. Matrix and vectors
	Rcpp::NumericVector yr(ParamVec);		// creates Rcpp vector from SEXP
	Rcpp::NumericMatrix Xr(Xs);			// creates Rcpp matrix from SEXP
	Rcpp::NumericMatrix G_f(FixedLinks);		// creates Rcpp matrix from SEXP
	int n = Xr.nrow(), k = Xr.ncol();
	arma::mat X_C(Xr.begin(), n, k, false);   	// reuses memory and avoids extra copy
	arma::colvec paramVec_C(yr.begin(), yr.size(), false);
	arma::mat Gamma_fixed(G_f.begin(), n, n, false);   	// reuses memory and avoids extra copy
// 	arma::umat uGamma_fixed = arma::conv_to<arma::umat>::from(Gamma_fixed);

	// .. String
	std::string ResFolder_C = Rcpp::as<std::string>(ResFolder);
	// .. ====================================================================================
// 	ProfilerStart("../../ProfileAR1");
	// .. Call mcmc function
	AR1_c(ResFolder_C, X_C, paramVec_C, Gamma_fixed);
// 	ProfilerStop();
	return 0;
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall

}



extern "C" SEXP callSplines(SEXP Xs, SEXP ResFolder, SEXP ParamVec, SEXP FixedLinks) {
  try {
	// .. Transform parameters to c format ===================================================	
	// .. Matrix and vectors
	Rcpp::NumericVector yr(ParamVec);		// creates Rcpp vector from SEXP
	Rcpp::NumericMatrix Xr(Xs);			// creates Rcpp matrix from SEXP
	Rcpp::NumericMatrix G_f(FixedLinks);		// creates Rcpp matrix from SEXP
	int n = Xr.nrow(), k = Xr.ncol();
	arma::mat X_C(Xr.begin(), n, k, false);   	// reuses memory and avoids extra copy
	arma::colvec paramVec_C(yr.begin(), yr.size(), false);
	arma::mat Gamma_fixed(G_f.begin(), n, n, false);   	// reuses memory and avoids extra copy

	// .. String
	std::string ResFolder_C = Rcpp::as<std::string>(ResFolder);
	// .. ====================================================================================
	
	// .. Call mcmc function
	PSplines_c(ResFolder_C, X_C, paramVec_C, Gamma_fixed);
	return 0;
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}


extern "C" SEXP callGauss_Error(SEXP Xs, SEXP ResFolder, SEXP ParamVec, SEXP FixedLinks) {
  try {
	// .. Transform parameters to c format ===================================================	
	// .. Matrix and vectors
	Rcpp::NumericVector yr(ParamVec);		// creates Rcpp vector from SEXP
	Rcpp::NumericMatrix Xr(Xs);			// creates Rcpp matrix from SEXP
	Rcpp::NumericMatrix G_f(FixedLinks);		// creates Rcpp matrix from SEXP
	int n = Xr.nrow(), k = Xr.ncol();
	arma::mat X_C(Xr.begin(), n, k, false);   	// reuses memory and avoids extra copy
	arma::colvec paramVec_C(yr.begin(), yr.size(), false);
	arma::mat Gamma_fixed(G_f.begin(), n, n, false);   	// reuses memory and avoids extra copy

	// .. String
	std::string ResFolder_C = Rcpp::as<std::string>(ResFolder);
	// .. ====================================================================================
	
	// .. Call mcmc function
	Error_Gauss_c(ResFolder_C, X_C, paramVec_C, Gamma_fixed);
	return 0;
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}



extern "C" SEXP callStudent_Error(SEXP Xs, SEXP ResFolder, SEXP ParamVec, SEXP FixedLinks)  {
  try {
	// .. Transform parameters to c format ===================================================	
	// .. Matrix and vectors
	Rcpp::NumericVector yr(ParamVec);		// creates Rcpp vector from SEXP
	Rcpp::NumericMatrix Xr(Xs);			// creates Rcpp matrix from SEXP
	Rcpp::NumericMatrix G_f(FixedLinks);		// creates Rcpp matrix from SEXP
	int n = Xr.nrow(), k = Xr.ncol();
	arma::mat X_C(Xr.begin(), n, k, false);   	// reuses memory and avoids extra copy
	arma::colvec paramVec_C(yr.begin(), yr.size(), false);
	arma::mat Gamma_fixed(G_f.begin(), n, n, false);   	// reuses memory and avoids extra copy

	// .. String
	std::string ResFolder_C = Rcpp::as<std::string>(ResFolder);
	// .. ====================================================================================
	// .. Call mcmc function
	Error_Student_c(ResFolder_C, X_C, paramVec_C, Gamma_fixed);
	return 0;
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}


extern "C" SEXP readLargeFileGetMean(SEXP Filename) {
  try {
	// .. Transform parameters to c format ===================================================	
	// .. Matrix and vectors
	arma::vec PostMeans;

	// .. String
	std::string Filename_C = Rcpp::as<std::string>(Filename);
	// .. ====================================================================================
	// .. Call mcmc function
	getPosteriorMeanFromFile(Filename_C, PostMeans);

	return Rcpp::wrap(PostMeans);
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}


extern "C" SEXP readGamma_getMean_numParents(SEXP Filename, SEXP FixedLinks) {
  try {
	// .. Transform parameters to c format ===================================================	
	// .. Matrix and vectors
	Rcpp::NumericMatrix G_f(FixedLinks);		// creates Rcpp matrix from SEXP
	int n = G_f.nrow();
	arma::mat Gamma_fixed(G_f.begin(), n, n, false);   	// reuses memory and avoids extra copy
	// .. String
	std::string Filename_C = Rcpp::as<std::string>(Filename);
	// .. ====================================================================================
	arma::mat numParents;
	arma::vec PostMeans;

	// .. Call mcmc function
        getPosteriorMeanFromFile_withNumParents(Filename_C, PostMeans, numParents, n, Gamma_fixed);	
	return Rcpp::List::create(PostMeans, numParents);
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}


