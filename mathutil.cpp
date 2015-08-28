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
#include "vector.hpp"
#include "logger.hpp"

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
// to mmin. This uses downward recursion, using the analytical relation of
// Mamedov, J. Math. Chem., 36, 3, 2004: "On the evaluation of Boys functions
// usign downward recursion relation". It is a very simple, efficient, and
// arbitrarily accurate (to within machine precision) method.
// DIGITS is the number of significant digits desired.
Vector boys(double x, int mmax, int mmin, int DIGITS)
{
  // Initialise return array
  Vector rvec(mmax - mmin + 1); // +1 as need to include both mmax and mmin
  double PRECISION = std::pow(10, -DIGITS);

  if(x < PRECISION) { // x ~ 0, so no need for recursion

    // x = 0 implies F_m(x) = 1/(2m+1) 
    for (int i = mmin; i <mmax + 1; i++){
      rvec[i-mmin] = 1.0/(2.0*(double)(i) + 1.0);
    }
    
  } else {

    // Declare starting value
    double f0 = std::sqrt(M_PI)/(2.0*sqrt(x));
  
    // Find the starting m for recursion from relation in Mamedov
    double dmmax = (double) mmax;
    double logval = std::log((fabs(dmmax - x) < PRECISION ? dmmax : dmmax/x));
    int mt = std::ceil( ((double)(DIGITS)/fabs(logval)) ) + mmax;

    // Come down the chain starting at mt, with F(mt) = f0
    // using standard boys downwards recursion relation:
    // F_m = 1/(2m+1)*(2xF_(m+1) + exp(-x))
    for (int i = mt; i > mmax; i--){
      f0 = ( 1.0/(2*(double)(i)) ) * ( 2*x*f0 + std::exp(-1.0*x) );
    }
    
    // Now carry on doing the same thing, but storing the results
    // in the return vector
    for (int i = mmax; i > mmin - 1; i--){
      f0 = ( 1.0/(2*(double)(i))) * ( 2*x*f0 + std::exp(-1.0*x) );
      rvec[i-mmin] = f0;
    }
  }
  return rvec;
}
  
