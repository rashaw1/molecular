/* 
 *    PURPOSE: To implement mathutil.hpp, a library of maths functions needed
 *             by various other programs.
 *
 *    DATE         AUTHOR          CHANGES
 *    =======================================================================
 *    27/08/15     Robert Shaw     Original Code
 *    28/08/15     Robert Shaw     Implemented Boys function.
 */

#include <cmath>
#include "mathutil.hpp"
#include "vector.hpp"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/constants/constants.hpp>

// Factorial and double factorial functions
long int fact(int i)
{
  long int rval = (i==0 ? 1 : i);
  for(int j = i-1; j > 0; j--){
    rval *= j;
  }
  return rval;
}

long int fact2(int i)
{
  long int rval = (i==0 ? 1 : i);
  for (int j = i-2; j > 0; j-=2){
    rval = rval*j;
  }
  return rval;
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
  Vector rvec(mmax - mmin + 1); // +1 as need to include both mmax and mmin

  if(x < PRECISION) { // x ~ 0, so no need for recursion

    // x = 0 implies F_m(x) = 1/(2m+1) 
    for (int i = mmin; i <mmax + 1; i++){
      rvec[i-mmin] = 1.0/(2.0*(double)(i) + 1.0);
    }
    
  } else {
    // For F_0, use erf
    if (mmin == 0){
      rvec[0] = std::sqrt(boost::math::pi()/(4.0*x))*boost::math::erf(std::sqrt(x));
    }

    // Calculate F_mmax value
    rvec[mmax-mmin] = (1.0/(2.0*std::pow(x, mmax+0.5))) * boost::math::tgamma_lower(mmax+0.5, x);
    
    // Now use downward recursion relation to calculate all others
    int m = (mmin == 0 ? 1 : mmin);
    for (int i = mmax-1; i > m; i--){
      rvec[i-mmin] = (1.0/(2.0*(double)(i) + 1.0)) * (2.0*x*rvec[i-mmin+1] + std::exp(-1.0*x));
    } 
  }
  return rvec;
}
  
