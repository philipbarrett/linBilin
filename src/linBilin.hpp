/**************************************************************************************
* linBilin.hpp                                                                        *
* Interface to linBilin.cpp                                                           *
* Philip Barrett, Chicago                                                             *
* Created: 25jul2015                                                                  *
*                                                                                     *
***************************************************************************************/

#ifndef HPP_LINBILIN
#define HPP_LINBILIN

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>  
#include <iostream> 

using namespace Rcpp ;
using namespace arma ;

double lin( NumericVector XX, NumericVector YY, double new_x, int length ) ;
    // Linear approximation
double approx_fast_arma( const arma::vec& XX, const arma::vec& YY, 
                          double new_x, int n_x ) ;
    // Armadillo interface
double approx_fast_arma_row( const arma::vec& XX, const arma::rowvec& YY, double new_x, int n_x ) ;
double approx_fast_arma_row2( const arma::rowvec& XX, const arma::vec& YY, double new_x, int n_x ) ;
double approx_fast_arma_row_row( const arma::rowvec& XX, const arma::rowvec& YY, 
                                  double new_x, int n_x ) ;
    // Row versions

double bilin_core( const arma::vec& x, const arma::vec& y, 
                          const arma::vec& vf, const arma::vec& pt ) ;
    // The heart of the bilinear interpolation algorithm
arma::umat bilin_idx_xy( const arma::vec& x, const arma::vec& y, 
                            const arma::vec& pt, int n_x, int n_y ) ;
arma::uvec bilin_idx_f( arma::umat m_idx, int n_x, int n_y ) ;
    // Index recovery prior to interpolation
double bilin( const arma::vec& pt, const arma::vec& x, const arma::vec& y, 
                  const arma::vec& f, int n_x, int n_y ) ;
    // Bilinear interpolation
double trilin( const arma::vec& pt, const arma::vec& x, const arma::vec& y, 
            const arma::vec& z, const arma::mat& m_f, int n_x, int n_y, int n_z ) ;
    // Trilinear interpolation
double trilin_z( const arma::vec& pt, const arma::vec& x, const arma::vec& y, 
                        const arma::mat& z, const arma::mat& m_f, 
                        int n_x, int n_y, int n_z ) ;
    // Trilinear interpolation with changing z
double trilin_inv( const arma::vec& pt, const arma::vec& x, const arma::vec& y, 
                      const arma::vec& z, const arma::mat& m_f, int n_x, int n_y, int n_z ) ;
    // Inverse trilinear interpolation
#endif