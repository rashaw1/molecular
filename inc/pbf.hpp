/*
 *
 *     PURPOSE: To define the class PBF representing a
 *              primitive cartesian gaussian basis function.
 *
 *     class PBF:
 *              data: exponent - the gaussian exponent of the basis function
 *                    norm - the normalisation constant
 *                    lx, ly, lz - the angular momentum quantum numbers in the 
 *                                 cartesian directions
 *              accessors: all have get routines
 *              routines:
 *                    normalise - calculate the normalisation constant
 *
 *     DATE        AUTHOR            CHANGES
 *     ====================================================================
 *     27/08/15    Robert Shaw       Original code.
 *
 */
 #ifndef PBFHEADERDEF
 #define PBFHEADERDEF
 
 class PBF
{
private:
  double exponent, norm;
  int lx, ly, lz;
public:
  // Constructors
  // Need to specify an exponent, e, and the angular momentum quantum numbers
  PBF() : exponent(0.0), norm(0.0), lx(0), ly(0), lz(0) { } // Default constructor
  PBF(double e, int l1, int l2, int l3);
  PBF(const PBF& other); // Copy constructor
  // Accessors
  double getExponent() const { return exponent; }
  double getNorm() const { return norm; }
  int getLnum() const { return lx+ly+lz; }
  int getLx() const { return lx; }
  int getLy() const { return ly; }
  int getLz() const { return lz; }
  // Routines
  void normalise();
  // Overloaded operators
  PBF& operator=(const PBF& other);
};
 
 #endif
