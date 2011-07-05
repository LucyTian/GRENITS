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
#include <iostream>
#include <cstdio>
#include <string>

extern "C" SEXP callAR1(SEXP Xs, SEXP ResFolder, SEXP ParamVec) {
  try {
	// .. Transform parameters to c format ===================================================	
	// .. Matrix and vectors
	Rcpp::NumericVector yr(ParamVec);		// creates Rcpp vector from SEXP
	Rcpp::NumericMatrix Xr(Xs);			// creates Rcpp matrix from SEXP
	int n = Xr.nrow(), k = Xr.ncol();
	arma::mat X_C(Xr.begin(), n, k, false);   	// reuses memory and avoids extra copy
	arma::colvec paramVec_C(yr.begin(), yr.size(), false);

	// .. String
	std::string ResFolder_C = Rcpp::as<std::string>(ResFolder);
	// .. ====================================================================================
	
	// .. Call mcmc function
	AR1_c(ResFolder_C, X_C, paramVec_C);
	return 0;
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall

}



extern "C" SEXP callSplines(SEXP Xs, SEXP ResFolder, SEXP ParamVec) {
  try {
	// .. Transform parameters to c format ===================================================	
	// .. Matrix and vectors
	Rcpp::NumericVector yr(ParamVec);		// creates Rcpp vector from SEXP
	Rcpp::NumericMatrix Xr(Xs);			// creates Rcpp matrix from SEXP
	int n = Xr.nrow(), k = Xr.ncol();
	arma::mat X_C(Xr.begin(), n, k, false);   	// reuses memory and avoids extra copy
	arma::colvec paramVec_C(yr.begin(), yr.size(), false);

	// .. String
	std::string ResFolder_C = Rcpp::as<std::string>(ResFolder);
	// .. ====================================================================================
	
	// .. Call mcmc function
	PSplines_c(ResFolder_C, X_C, paramVec_C);
	return 0;
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}


extern "C" SEXP callGauss_Error(SEXP Xs, SEXP ResFolder, SEXP ParamVec) {
  try {
	// .. Transform parameters to c format ===================================================	
	// .. Matrix and vectors
	Rcpp::NumericVector yr(ParamVec);		// creates Rcpp vector from SEXP
	Rcpp::NumericMatrix Xr(Xs);			// creates Rcpp matrix from SEXP
	int n = Xr.nrow(), k = Xr.ncol();
	arma::mat X_C(Xr.begin(), n, k, false);   	// reuses memory and avoids extra copy
	arma::colvec paramVec_C(yr.begin(), yr.size(), false);

	// .. String
	std::string ResFolder_C = Rcpp::as<std::string>(ResFolder);
	// .. ====================================================================================
	
	// .. Call mcmc function
	Error_Gauss_c(ResFolder_C, X_C, paramVec_C);
	return 0;
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}



extern "C" SEXP callStudent_Error(SEXP Xs, SEXP ResFolder, SEXP ParamVec) {
  try {
	// .. Transform parameters to c format ===================================================	
	// .. Matrix and vectors
	Rcpp::NumericVector yr(ParamVec);		// creates Rcpp vector from SEXP
	Rcpp::NumericMatrix Xr(Xs);			// creates Rcpp matrix from SEXP
	int n = Xr.nrow(), k = Xr.ncol();
	arma::mat X_C(Xr.begin(), n, k, false);   	// reuses memory and avoids extra copy
	arma::colvec paramVec_C(yr.begin(), yr.size(), false);

	// .. String
	std::string ResFolder_C = Rcpp::as<std::string>(ResFolder);
	// .. ====================================================================================
	// .. Call mcmc function
	Error_Student_c(ResFolder_C, X_C, paramVec_C);
	return 0;
    } catch( std::exception &ex ) {
	forward_exception_to_r( ex );
    } catch(...) { 
	::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}
