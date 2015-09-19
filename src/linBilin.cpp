/**************************************************************************************
* linBilin.cpp                                                                        *
* Code to do fast linear and bilinear approximation                                   *
* Philip Barrett, Chicago                                                             *
* Created: 25jul2015                                                                  *
*                                                                                     *
***************************************************************************************/

#include "linBilin.hpp"

/***********************************************************************/
/** 1. LINEAR APPROXIMATION                                           **/
/***********************************************************************/

// [[Rcpp::export]]
double lin( NumericVector x, NumericVector y, double pt, int n_x ){
// Linear approximation of a vector-function
  if( pt < x[0] ){
    return( y[0] ) ;
  }
  if( pt >= x[n_x-1] ){
    return( y[n_x-1] ) ;
  }
  
  int idx = 0 ;
      // Initialize index
  while ( x[idx] <= pt && idx < n_x ){
    idx++ ;
  } // Find the location of pt in x

  double theta = ( pt - x[idx-1] ) / ( x[idx] - x[idx-1] ) ;
      // Fraction of distance between two x-values
  return theta * y[idx] + ( 1 - theta ) * y[idx-1] ;
}

double approx_fast_arma( const arma::vec& XX, const arma::vec& YY, double new_x, int n_x ){
// Linear approximation of a vector-function
  
  if( new_x <= XX[0] ){
    return YY[0] ;
  }
  if ( new_x >= XX[n_x-1] ){
    return YY[n_x-1] ;
  }
  
  uvec v_idx = find( XX > new_x, 1, "first" ) ;
  int idx = v_idx[0] ;
      // Index
  double out = ( ( new_x - XX[idx-1] ) * YY[idx] + ( XX[idx] - new_x ) * YY[idx-1] ) / ( XX[idx] - XX[idx-1] ) ;
      // Fraction of distance between two x-values
  return out ;
}

double approx_fast_arma_row( const arma::vec& XX, const arma::rowvec& YY, double new_x, int n_x ){
// Linear approximation of a vector-function
  
  if( new_x <= XX[0] ){
    return YY[0] ;
  }
  if ( new_x >= XX[n_x-1] ){
    return YY[n_x-1] ;
  }
  
  uvec v_idx = find( XX > new_x, 1, "first" ) ;
  int idx = v_idx[0] ;
      // Index
  double out = ( ( new_x - XX[idx-1] ) * YY[idx] + ( XX[idx] - new_x ) * YY[idx-1] ) / ( XX[idx] - XX[idx-1] ) ;
      // Fraction of distance between two x-values
  return out ;
}

double approx_fast_arma_row_row( const arma::rowvec& XX, const arma::rowvec& YY, 
                                  double new_x, int n_x ){
// Linear approximation of a vector-function
  
  if( new_x <= XX[0] ){
    return YY[0] ;
  }
  if ( new_x >= XX[n_x-1] ){
    return YY[n_x-1] ;
  }
  
  uvec v_idx = find( XX > new_x, 1, "first" ) ;
  int idx = v_idx[0] ;
      // Index
  double out = ( ( new_x - XX[idx-1] ) * YY[idx] + ( XX[idx] - new_x ) * YY[idx-1] ) / ( XX[idx] - XX[idx-1] ) ;
      // Fraction of distance between two x-values
  return out ;
}

double approx_fast_arma_row2( const arma::rowvec& XX, const arma::vec& YY, double new_x, int n_x ){
// Linear approximation of a vector-function
  
  if( new_x <= XX[0] ){
    return YY[0] ;
  }
  if ( new_x >= XX[n_x-1] ){
    return YY[n_x-1] ;
  }
  
  uvec v_idx = find( XX > new_x, 1, "first" ) ;
  int idx = v_idx[0] ;
      // Index
  double out = ( ( new_x - XX[idx-1] ) * YY[idx] + ( XX[idx] - new_x ) * YY[idx-1] ) / ( XX[idx] - XX[idx-1] ) ;
      // Fraction of distance between two x-values
  return out ;
}

/***********************************************************************/
/** 2. BILINEAR APPROXIMATION                                         **/
/***********************************************************************/

double bilin_core( const arma::vec& x, const arma::vec& y, 
                          const arma::vec& vf, const arma::vec& pt ){
// Bilinear approximation of a vector-function at a point pt in a square with x 
// in (x1,x2), y in (y1,y2), and vf= ( f(x1,y1), f(x1,y2), f(x2,y1), f(x2,y2) ).
// See:
// https://en.wikipedia.org/wiki/Bilinear_interpolation#Alternative_algorithm 
// for further details
  
  // Check that the point is not exact
  if( x[0] == pt[0] ){
    if( y[0] == pt[1] ){
      return vf[0] ;
    }
    else if( y[1] == pt[1] ){
      return vf[1] ;
    }
  }
  else if( x[1] == pt[0] ){
    if( y[0] == pt[1] ){
      return vf[2] ;
    }
    else if( y[1] == pt[1] ){
      return vf[3] ;
    }
  }
  
  // Check that the point is not on the boundary and do 1D interpolation
  if( x[0] == x[1] ){
    if( y[0] == y[1] ){
      return( vf[0] ) ;
    }
    return ( vf[0] * ( y[1] - pt[1] ) + vf[1]  * ( pt[1] - y[0] ) ) / ( y[1] - y[0] ) ;
  }
  if( y[0] == y[1] ){
    return ( vf[0] * ( x[1] - pt[0] ) + vf[2]  * ( pt[0] - x[0] ) ) / ( x[1] - x[0] ) ;
  }
  
  
  // Solve on the interior
  mat premult( 4, 4 ) ;
  premult << 1 << x[0] << y[0] << x[0] * y[0] << endr
          << 1 << x[0] << y[1] << x[0] * y[1] << endr
          << 1 << x[1] << y[0] << x[1] * y[0] << endr
          << 1 << x[1] << y[1] << x[1] * y[1] << endr ;
      // The premultiplying matrix in the solution for the approximation coefficients
  
//  Rcout << "premult: \n" << premult << std::endl ;
  
  vec coefs = solve( premult, vf ) ;
      // Coefficients of the approximation
  return coefs[0] + coefs[1] * pt[0] + coefs[2] * pt[1] + coefs[3] * pt[0] * pt[1] ;
}


arma::umat bilin_idx_xy( const arma::vec& x, const arma::vec& y, 
                          const arma::vec& pt, int n_x, int n_y ){
// Given a point pt, identifies the indices of the four points in the x*y grid 
// which contain pt.  The return is the index in the table that lists the grid
// points first by x and then by y.

  umat idx(2,2) ;
      // The matrix of indices in x & y
  
  int i_x_low_idx, i_x_high_idx ;
  if( pt[0] <= x[0] ){
    i_x_low_idx = 0 ;
    i_x_high_idx = 0 ;
  }
  else if( pt[0] > x[n_x-1] ){
    i_x_low_idx = n_x - 1 ;
    i_x_high_idx = n_x - 1 ;
  }
  else{
    uvec v_x_low_idx = find( x < pt[0], 1, "last" ) ;
    i_x_low_idx = v_x_low_idx(0) ;
    i_x_high_idx = v_x_low_idx(0) + 1 ;
  } // The x index
  
  int i_y_low_idx, i_y_high_idx ;
  if( pt[1] <= y[0] ){
    i_y_low_idx = 0 ;
    i_y_high_idx = 0 ;
  }
  else if( pt[1] > y[n_y-1] ){
    i_y_low_idx = n_y - 1 ;
    i_y_high_idx = n_y - 1 ;
  }
  else{
    uvec v_y_low_idx = find( y < pt[1], 1, "last" ) ;
    i_y_low_idx = v_y_low_idx(0) ;
    i_y_high_idx = v_y_low_idx(0) + 1 ;
  } // The y index

  idx << i_x_low_idx << i_y_low_idx << endr
      << i_x_high_idx << i_y_high_idx ;
  
  return idx ;
      // NB: Returns c-style (zero-based) indices
}

arma::uvec bilin_idx_f( arma::umat m_idx, int n_x, int n_y ){
// Returns the indices of f in the table ordered first by x then by y

  uvec idx(4) ;
  idx[0] = m_idx(0,0) * n_y + m_idx(0,1) ;
  idx[1] = m_idx(0,0) * n_y + m_idx(1,1) ;
  idx[2] = m_idx(1,0) * n_y + m_idx(0,1) ;
  idx[3] = m_idx(1,0) * n_y + m_idx(1,1) ;
      // Work out the locations in the table
      
  return idx ;
}

// [[Rcpp::export]]
double bilin( const arma::vec& pt, const arma::vec& x, const arma::vec& y, 
                  const arma::vec& f, int n_x, int n_y ){
// Finds the bilinear interpolation of the values f over the grid x * y, where f
// is listed by x first and y second
  
  umat idx_xy = bilin_idx_xy( x, y, pt, n_x, n_y ) ;
  uvec idx_f = bilin_idx_f( idx_xy, n_x, n_y ) ;
      // The indices of the surrounding box and values in f
  
  vec x_bilin = x.elem( idx_xy.col( 0 ) ) ;
  vec y_bilin = y.elem( idx_xy.col( 1 ) ) ;
      // The values of x and y to pass to the bilinear interpolator
  vec f_bilin = f( idx_f ) ;
      // The corresponding values of f
  return bilin_core( x_bilin, y_bilin, f_bilin, pt ) ;
}


/***********************************************************************/
/** 3. TRILINEAR APPROXIMATION                                        **/
/***********************************************************************/

// [[Rcpp::export]]
double trilin( const arma::vec& pt, const arma::vec& x, const arma::vec& y, 
               const arma::vec& z, const arma::mat& m_f, int n_x, int n_y, int n_z ){
// Performs trilinear interpolation, estimating the value of f at the point pt 
// on the grid x*y*z, where m_f is a vector of the function evaluated at x then
// y (vertically) and z (horizontally).
  
  vec pt_xy ;
  pt_xy << pt(0) << pt(1) << endr ;
  double d_z = pt(2) ;
  
  umat idx_xy = bilin_idx_xy( x, y, pt_xy, n_x, n_y ) ;
  uvec idx_f = bilin_idx_f( idx_xy, n_x, n_y ) ;
      // The indices of the surrounding box and values in f
      
  vec vf(4) ;
  for( int i = 0; i < 4 ; i++ ){
    rowvec vv = m_f.row( idx_f(i) ) ;
    vf(i) = approx_fast_arma_row( z, m_f.row( idx_f(i) ), d_z, n_z ) ;
  }  // Linear approx of the function in the z dimension

  vec x_bilin = x.elem( idx_xy.col( 0 ) ) ;
  vec y_bilin = y.elem( idx_xy.col( 1 ) ) ;
      // The values of x and y to pass to the bilinear interpolator

  return bilin_core( x_bilin, y_bilin, vf, pt_xy ) ;
}

// [[Rcpp::export]]
double trilin_z( const arma::vec& pt, const arma::vec& x, const arma::vec& y, 
                        const arma::mat& z, const arma::mat& m_f, 
                        int n_x, int n_y, int n_z ){
// Performs trilinear interpolation, when the grid of points in the z-dimension
// is not a single fixed vector, but rather a matrix of changing values
  
  vec pt_xy ;
  pt_xy << pt(0) << pt(1) << endr ;
  double d_z = pt(2) ;
  
  umat idx_xy = bilin_idx_xy( x, y, pt_xy, n_x, n_y ) ;
  uvec idx_f = bilin_idx_f( idx_xy, n_x, n_y ) ;
      // The indices of the surrounding box and values in f
      
  vec vf(4) ;
  for( int i = 0; i < 4 ; i++ ){
    vf(i) = approx_fast_arma_row_row( z.row( idx_f(i) ), m_f.row( idx_f(i) ), d_z, n_z ) ;
  }  // Linear approx of the function in the z dimension

  vec x_bilin = x.elem( idx_xy.col( 0 ) ) ;
  vec y_bilin = y.elem( idx_xy.col( 1 ) ) ;
      // The values of x and y to pass to the bilinear interpolator

  return bilin_core( x_bilin, y_bilin, vf, pt_xy ) ;
}


// [[Rcpp::export]]
double trilin_inv( const arma::vec& pt, const arma::vec& x, const arma::vec& y, 
                      const arma::vec& z, const arma::mat& m_f, int n_x, int n_y, int n_z ){
// Performs invesre trilinear interpolation, estimating the value of z at the 
// point on the grid x*y*z such that f(x,y,z)=q, where q is the third entry of 
// pt, and m_f is a vector of the function evaluated at x then y (vertically)
// and z (horizontally).
  
  vec pt_xy(2) ;
  pt_xy << pt(0) << pt(1) << endr ;
  double d_q = pt(2) ;
  
  umat idx_xy = bilin_idx_xy( x, y, pt_xy, n_x, n_y ) ;
  uvec idx_f = bilin_idx_f( idx_xy, n_x, n_y ) ;
      // The indices of the surrounding box and values in f
      
  vec v_z(4) ;
  for( int i = 0; i < 4 ; i++ ){
    v_z(i) = approx_fast_arma_row2( m_f.row( idx_f(i) ), z, d_q, n_z ) ;
  }  // Linear approx of the function in the z dimension
  
  vec x_bilin = x.elem( idx_xy.col( 0 ) ) ;
  vec y_bilin = y.elem( idx_xy.col( 1 ) ) ;
      // The values of x and y to pass to the bilinear interpolator
  
  double z_new = bilin_core( x_bilin, y_bilin, v_z, pt_xy ) ;
  double z_old, z_1, z_2, f_1, f_2, diff ;
  uvec v_z_1_idx, v_z_2_idx ;
  int i_z_1_idx, i_z_2_idx ;
  vec vf_1(4), vf_2(4), test_pt(3) ;
  bool flat = false ;
  
//  Rcout << "z_new = " << z_new << std::endl ;
  
  int i = 0 ;
  int maxit = 10 ;
  double tol = 1e-08 ;
  diff = 2 * tol ;
      // Initialize variables for iterative search
  
  while( ( i < maxit ) & ( diff > tol ) ){
    z_old=z_new ;
        // Update the guess
    v_z_1_idx = find( z < z_old, 1, "last" ) ;
    v_z_2_idx = v_z_1_idx + 1 ;
        // Find the indices of the current z
    i_z_1_idx = ( z_old <= z(0) ? 0 : v_z_1_idx(0) ) ;
    i_z_2_idx = ( z_old > z(n_z-1) ? n_z-1 : v_z_2_idx(0) ) ;
        // Convert to integers
        
//    Rcout << "i_z_1_idx = " << i_z_1_idx << std::endl ;
//    Rcout << "i_z_2_idx = " << i_z_2_idx << std::endl ;
        
    for( int j = 0 ; j < 4 ; j++ ){
      vf_1( j ) = m_f( idx_f( j ), i_z_1_idx ) ;
      vf_2( j ) = m_f( idx_f( j ), i_z_2_idx ) ;
    } // Find the value of f on the boundaries of the box containing z_old
    
    f_1 = bilin_core( x_bilin, y_bilin, vf_1, pt_xy ) ;
    f_2 = bilin_core( x_bilin, y_bilin, vf_2, pt_xy ) ;
        // Compute f_1=f(x,y,z_1) and f_2=f(x,y,z_2)
        
    z_1 = z( i_z_1_idx ) ;
    z_2 = z( i_z_2_idx ) ;
        // Update the values of z
        
    if( f_1 == f_2 ){
      z_new = .5 * ( z_1 + z_2 ) ;
      flat = true ;
          // If we have a flat spot, just pick the point in the middle.
    }
    else{
      z_new = ( ( d_q - f_1 ) * z_2 + ( f_2 - d_q ) * z_1 ) / ( f_2 - f_1 ) ;
          // Solves:
          // ( z - z_1 ) / (z_2-z_1) * f_2 + (z_2-z) / (z_2-z_1) * f_1 = q
      
    }
    test_pt << pt(0) << pt(1) << z_new << endr ;
    diff = fabs( d_q - trilin( test_pt, x, y, z, m_f, n_x, n_y, n_z ) ) ;
    i++ ;
      // Update i.  Should usually only need one iteration to solve.
  }
  
  if( !flat ){
    if( diff > tol ){
      Rcout << "WARNING: TRILINEAR INTERPOLIATION INVERSION FAILING.  ERROR = " << 
        round( diff * 10000.0 ) / 10000.0 << std::endl ;
    }
  }
  
return z_new ;
}