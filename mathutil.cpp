/* 
 *    PURPOSE: To implement mathutil.hpp, a library of maths functions needed by 
 *             various other programs.
 *
 *    DATE         AUTHOR          CHANGES
 *    ==========================================================================
 *    27/08/15     Robert Shaw     Original Code
 *
 */

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
