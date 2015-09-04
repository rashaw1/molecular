/*
 *    PURPOSE: To implement bf.hpp, defining class BF, a contracted
 *             cartesian gaussian basis function.
 *   
 *    DATE            AUTHOR              CHANGES
 *    =====================================================================
 *    27/08/15        Robert Shaw         Original code.
 *
 */

#include "bf.hpp"

BF::BF(Vector& c, int l1, int l2, int l3, Vector& exps, Vector& indices)
{
  coeffs = c;
  ids = indices;
  lx = l1; ly = l2; lz = l3;
  
  // Construct the set of pbfs
  int npbfs = exps.size();
  if (npbfs > 0){
    pbfs = new PBF[npbfs];
    for (int i = 0; i < npbfs; i++){
      PBF temp(exps(i), l1, l2, l3);
      pbfs[i] = temp;
    }
  } else {
    pbfs = NULL;
  }
}

// Copy constructor
BF::BF(const BF& other)
{
  coeffs = other.coeffs;
  ids = other.ids;
  norm = other.norm;
  lx = other.lx; ly = other.ly; lz = other.lz;
  
  int npbfs = coeffs.size();
  if(npbfs > 0){
    pbfs = new PBF[npbfs];
    for (int i = 0; i < npbfs; i++){
      pbfs[i] = other.pbfs[i];
    }
  } else {
    pbfs = NULL;
  }
}

// Destructor
BF::~BF()
{
  if (coeffs.size() > 0){
    delete[] pbfs;
  }
}

// Accessors
Vector BF::getExps() const
{
  // Construct a vector of exps
  Vector exps(coeffs.size());
  for (int i = 0; i < coeffs.size(); i++){
    exps[i] = pbfs[i].getExponent();
  }
  return exps;
}

BF& BF::operator=(const BF& other)
{
  // If PBFs already exist, deallocate memory
  if(coeffs.size() > 0){
    delete[] pbfs;
  }
  
  // Assign attributes
  coeffs = other.coeffs;
  ids = other.ids;
  norm = other.norm;
  lx = other.lx; ly = other.ly; lz = other.lz;

  // Copy across PBFs
  int npbfs = coeffs.size();
  if (npbfs > 0){
    pbfs = new PBF[npbfs];
    for (int i = 0; i < npbfs; i++){
      pbfs[i] = other.pbfs[i];
    }
  }
  
  return *this;
}
