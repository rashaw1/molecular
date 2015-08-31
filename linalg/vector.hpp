/*
 *     PURPOSE: defines and implements class Vector, a vector of doubles 
 *
 * 
 *     DATE             AUTHOR                CHANGES
 *     =====================================================================
 *     13/08/15         Robert Shaw           Original code
 *     14/08/15         Robert Shaw           Added error throwing.
 *     20/08/15         Robert Shaw           Added outer product, angle,
 *                                            and sorting.
 *     26/08/15         Robert Shaw           Added triple and cross products.
 */

#ifndef VECTORHEADERDEF
#define VECTORHEADERDEF

#include "matrix.hpp"
#include "error.hpp"

class Vector
{
private:
  int n; // The number of elements
  double* v; // The elements themselves
  void cleanUp(); // Deallocates memory
public:
  // Constructors and destructor
  Vector() : n(0), v(NULL) {} // Default constructor, zero length vector
  Vector(int length); // Empty vector of length length
  Vector(int length, const double& a); // Vector with 'length' values, all a
  Vector(int length, const double* a); // Initialise vector to array a
  Vector(const Vector& u); // Copy constructor
  ~Vector(); // Destructor
  // Accessors
  int size() const { return n; } // Returns size of vector, n
  // Shaping functions
  void resize(int length); // Resizes the vector to length 'length',
                           // without preserving values
  void resizeCopy(int length); // Resizes and preserves values up to length
  void assign(int length, const double& a); // Resizes and sets elements to a
  void swap(int i, int j); // Swaps elements i and j
  // Overloaded operators
  double& operator[](int i); // Access value at index i
  double operator[](int i) const; // Return by value
  double operator()(int i) const; // Also return by value
  Vector& operator=(const Vector& u); // Set this = u
  // Unary operators
  Vector operator+() const; 
  Vector operator-() const;
  // Binary operators
  Vector operator+(const Vector& u) const;
  Vector operator-(const Vector& u) const;
  Vector& operator*=(const double& scalar) { return *this; } // Scalar multiplication
  Vector& operator*=(const Matrix& mat) { return *this; } // Vector x matrix
  // Intrinsic functions
  void print(double PRECISION = 1e-12) const; // Pretty prints the vector to primary ostream
  Vector sorted() const; // Returns a sorted copy of the vector
  void sort(); // Sorts the vector into ascending order, uses quicksort 
  // Friend functions
  friend double pnorm(const Vector& u, int p); // Returns the p-norm of u
  // Calculate the inner (dot) product of two vectors
  friend double inner(const Vector& u, const Vector& w);
  // Calculate the outer product of two vectors
  friend Matrix outer(const Vector& u, const Vector& w);
  // Return angle (in radians) between two vectors
  friend double angle(const Vector& u, const Vector& w);

  // !!!ONLY FOR 3D VECTORS!!!
  // Calculate the cross product of two vectors
  friend Vector cross(const Vector& u, const Vector& w);
  // Calculate the triple product of three vectors
  friend double triple(const Vector& u, const Vector& w, const Vector& z);
};

// Scalar multiplication

inline Vector operator*(const double& scalar, const Vector& v)
{
  int n = v.size();
  Vector rvec(n);
  for(int i = 0; i < n; i++){
    rvec[i] = v(i)*scalar;
  }
  return rvec;
}

inline Vector operator*(const Vector& v, const double& scalar)
{
  int n = v.size();
  Vector rvec(n);
  for (int i = 0; i < n; i++){
    rvec[i] = scalar*v(i);
  }
  return rvec;
}

// Vector x matrix and matrix x vector- will throw an error if wrong shapes

inline Vector operator*(const Vector& v, const Matrix& mat)
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


inline Vector operator*(const Matrix& mat, const Vector& v)
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


#endif
