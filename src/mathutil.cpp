/* 
 *    PURPOSE: To implement mathutil.hpp, a library of maths functions needed
 *             by various other programs.
 *
 *    DATE         AUTHOR          CHANGES
 *    =======================================================================
 *    27/08/15     Robert Shaw     Original Code
 *    28/08/15     Robert Shaw     Implemented Boys function.
 *    02/09/15     Robert Shaw     Implemented binom function.
 *    05/09/15     Robert Shaw     Implemented clebsch, sphernorm.
 *    07/09/15     Robert Shaw     Implemented formTransMat.
 *    12/09/15     Robert Shaw     Move vecxmatrix multiplication here.
 */

#include <cmath>
#include "mathutil.hpp"
#include "error.hpp"
#include "mvector.hpp"
#include "matrix.hpp"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/constants/constants.hpp>

// Factorial and double factorial functions
unsigned long int fact(int i)
{
  unsigned long int rval = (i<1 ? 1 : i);
  for(int j = i-1; j > 0; j--){
    rval *= j;
  }
  return rval;
}

unsigned long int fact2(int i)
{
  unsigned long int rval = (i<1 ? 1 : i);
  for (int j = i-2; j > 0; j-=2){
    rval = rval*j;
  }
  return rval;
}

void factArray(int i, double *values) {
  values[0] = 1.0;
  if ( i > 0 ) {
    values[1] = 1.0;
    for (int j = 2; j <= i; j++) values[j] = values[j-1]*j;
  }				
}

void fact2Array(int i, double *values) {
  values[0] = 1.0;
  if ( i > 0 ) {
    values[1] = 1.0;
    for (int j = 2; j <= i; j++) values[j] = values[j-2]*j;
  }
}

// Calculate the boys function F_m(x) for a range of values of m from mmax
// to mmin.
// First, the mmax value is calculated using the relation to the (full)
// incomplete gamma function, then a downward recursion is used for all
// values of m down to mmin. If F_0 is needed, this is calculated by a
// call to the error function instead (speedier). 
// These relations can all be found in Helgaker, Jorgensen, Olsen,
// "Molecular Electronic Structure Theory", Chapter 9, section 8, 
// pages 366-371, (in the 2012 paperback edition).
// Using boost libraries for their high accuracy and efficiency.

Vector boys(double x, int mmax, int mmin, double PRECISION)
{
  // Initialise return array
  Vector rvec(mmax - mmin + 1, 0.0); // +1 as need to include both mmax and mmin

  if(x < PRECISION) { // x ~ 0, so no need for recursion

    // x = 0 implies F_m(x) = 1/(2m+1) 
    for (int i = mmin; i <mmax + 1; i++){
      rvec[i-mmin] = 1.0/(2.0*(double)(i) + 1.0);
    }
    
  } else {
    // For F_0, use erf
    if (mmin == 0){
      rvec[0] = std::sqrt(boost::math::constants::pi<double>()/(4.0*x))*boost::math::erf(std::sqrt(x));
    }

    // Calculate F_mmax value
    rvec[mmax-mmin] = (1.0/(2.0*std::pow(x, mmax+0.5))) * boost::math::tgamma_lower(mmax+0.5, x);
    
    // Now use downward recursion relation to calculate all others
    int m = (mmin == 0 ? 1 : mmin);
    for (int i = mmax-1; i > m-1; i--){
      rvec[i-mmin] = (1.0/(2.0*(double)(i) + 1.0)) * (2.0*x*rvec[i-mmin+1] + std::exp(-1.0*x));
    } 
  }
  return rvec;
}

// Calculate the binomial coefficient (n m)T
// given by the formula:
// n! / [(n-m)!*m!]
unsigned long int binom(int n, int m)
{  
  unsigned long int rval;
  if (m >= n || m <= 0 || n <= 0){
    rval = 1;
  } else {
    unsigned long int numerator = n;
    // Set m to be the lesser of m, n-m
    m = (m < n-m ? m : n-m);
    if ( m == 1 ) {
      rval = n;
    } else {
      unsigned long int denominator = m;
      for (int i = m-1; i > 0; i--){
	numerator *= (n-i);
	denominator *= i;
      }
      rval =  numerator/denominator;
    }
  }
  return rval;
}

// Calculate the Clebsch-Gordan Coefficient C^{lm}_{tuv}
// according to the formula in Helgaker, Jorgensen, Olsen,
// Molecular Electronic Structure Theory, chapter 9 pg 338
double clebsch(int l, int m, int t, int u, double v)
{
  double vm = (m < 0 ? 0.5 : 0.0);
  int sign = ( ((int)(t+v-vm))%2 == 0 ? 1 : -1);
  double cval = sign*std::pow(0.25, t);

  // Get the binomial factors
  int mabs = std::abs(m);
  double ltb = binom(l, t);
  double lmtb = binom(l-t, mabs + t);
  double tub = binom(t, u);
  cval = cval*ltb*lmtb*tub*binom(mabs, (int)(2*v));
  return cval;
}

// Calculate the spherical normalisation Nlm
double sphernorm(int l, int m)
{
  int mabs = std::abs(m);
  double nval = std::pow(0.5, mabs)/((double)(fact(l)));
  
  // Get factorials
  double lpmfact = fact(l + mabs);
  double lmmfact = fact(l - mabs);
  double zerom = (m==0 ? 2.0 : 1.0);
  nval = nval*std::sqrt(2*lpmfact*lmmfact/zerom);
  return nval;
}

// Fill in a portion of the cart-to-spher. transformation matrix
// corresponding to quantum numbers l, m. Assumes canonical order.
void formTransMat(Matrix& mat, int row, int col, int l, int m){
  double temp; 
  switch(l){
  case 1: { // p-type, one of 3 identity 3-vectors
    switch(m){
    case -1: { mat(row, col+1) = 1.0; break; } // py
    case 0: { mat(row, col+2) = 1.0; break; }  // pz
    default: { mat(row, col) = 1.0;} // px
    }
    break;
  }
  case 2: { // d-type
    switch(m) {
    case -2: { mat(row, col+4) = 1.0; break; } // dxy
    case -1: { mat(row, col+1) = 1.0; break; } // dyz
    case 0: { // dz2 
      mat(row, col) = 1.0; 
      mat(row, col+3) = mat(row, col+5) = -0.5;
      break;
    }
    case 1: { mat(row, col+2) = 1.0; break; } // dxz
    default: { // d(x^2-y^2)
      temp = std::sqrt(3.0/4.0);
      mat(row, col+5) = temp;
      mat(row, col+3) = -temp;
    }
    }
    break;
  }
  case 3: { // f-type
    switch(m){
    case -3: { 
      mat(row, col+6) = std::sqrt(10)/4.0;
      mat(row, col+8) = -3.0*std::sqrt(2)/4.0;
      break;
    }
    case -2: { mat(row, col+4) = 1.0; break; }
    case -1: { 
      temp = std::sqrt(6.0/5.0);
      mat(row, col+1) = temp;
      mat(row, col+6) = -std::sqrt(6)/4.0;
      mat(row, col+8) = -temp/4.0;
      break;
    }
    case 0: {
      mat(row, col) = mat(row, col+3) = 1.0;
      mat(row, col+5) = -3.0/(2.0*std::sqrt(5));
      break;
    }
    case 1: {
      temp = std::sqrt(6.0/5.0);
      mat(row, col+2) = temp;
      mat(row, col+9) = -std::sqrt(6)/4.0;
      mat(row, col+7) = -temp/4.0;
      break;
    }
    case 2: {
      temp = std::sqrt(3.0/4.0);
      mat(row, col+5) = temp;
      mat(row, col+3) = -temp;
      break;
    }
    default: {
      mat(row, col+9) = std::sqrt(10)/4.0;
      mat(row, col+7) = -3.0*std::sqrt(2)/4.0;
    }
    }
    break;
  }
  case 4: { // g-type
    switch(m) {
    case -4:{
      temp = std::sqrt(5.0/4.0);
      mat(row, col+13) = temp;
      mat(row, col+11) = -temp;
      break;
    }
    case -3:{
      mat(row, col+6) = -1.0*std::sqrt(10.0)/4.0;
      mat(row, col+8) = 0.75*std::sqrt(2.0);
      break;
    }
    case -2:{
      temp = -1.0*std::sqrt(5.0/7.0)/2.0;
      mat(row, col+4) = 3.0*std::sqrt(1.0/7.0);
      mat(row, col+13) = temp;
      mat(row, col+11) = temp;
      break;
    }
    case -1:{
      temp = std::sqrt(2.0/7.0);
      mat(row, col+1) = temp*std::sqrt(5.0);
      mat(row, col+6) = -0.75*temp*std::sqrt(5.0);
      mat(row, col+8) = -0.75*temp;
      break;
    }
    case 0:{
      temp = -3.0*std::sqrt(3.0/35/0);
      mat(row, col) = 1.0;
      mat(row, col+14) = mat(row, col+10) = 3.0/8.0;
      mat(row, col+5) = mat(row, col+3) = temp;
      mat(row, col+12) = -0.25*temp;
      break;
    }
    case 1:{ 
      temp = std::sqrt(2.0/7.0);
      mat(row, col+2) = std::sqrt(5.0)*temp;
      mat(row, col+9) = -0.75*std::sqrt(5.0)*temp;
      mat(row, col+7) = -0.75*temp;
      break;
    }
    case 2:{
      temp = 1.5*std::sqrt(3.0/7.0);
      mat(row, col+5) = temp;
      mat(row, col+3) = -temp;
      temp = 0.25*std::sqrt(5.0);
      mat(row, col+14) = -temp;
      mat(row, col+10) = temp;
      break;
    }
    case 3:{
      temp = 0.25*std::sqrt(2.0);
      mat(row, col+9) = temp*std::sqrt(5.0);
      mat(row, col+7) = -3.0*temp;
      break;
    }
    default:{
      temp = std::sqrt(35.0)/8.0;
      mat(row, col+14) = mat(row, col+10) = temp;
      mat(row, col+12) = -0.75*std::sqrt(3.0);
    }     
    }
    break;
  }
  default: { // s-type
    mat(row, col) = 1.0;
  }
  }
}


// Vector x matrix and matrix x vector- will throw an error if wrong shapes                                                                                         
// Left and right vector x matrix multiplication functions

Vector lmultiply(const Vector& v, const Matrix& mat)
{
  // Assume left multiplication implies transpose                                                                                                                   
  int n = v.size();
  int cols = mat.ncols();
  int rows = mat.nrows();
  Vector rVec(cols); // Return vector should have dimension cols                                                                                                    
  // For this to work we require n = rows                                                                                                                           
  if (n != rows) {
    throw(Error("VECMATMULT", "Vector and matrix are wrong sizes to multiply."));
  } else { // Do the multiplication                                                                                                                                 
    for (int i = 0; i < cols; i++){
      Vector temp(rows); // Get each column of the matrix                                                                                                           
      temp = mat.colAsVector(i);
      rVec[i] = inner(v, temp);
    }
  }
  return rVec;
}


Vector rmultiply(const Matrix& mat, const Vector& v)
{
  int n = v.size();
  int cols = mat.ncols();
  int rows = mat.nrows();
  Vector rVec(rows); // Return vector should have dimension rows                                                                                                    
  // For this to work we require n = cols                                                                                                                           
  if (n != cols) {
    throw(Error("MATVECMULT", "Vector and matrix are wrong sizes to multiply."));
  } else { // Do the multiplication                                                                                                                                 
    for (int i = 0; i < rows; i++){
      Vector temp(cols);
      temp = mat.rowAsVector(i);
      rVec[i] = inner(v, temp);
    }
  }
  return rVec;
}
      
