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

class Matrix;

#include "error.hpp"
#include "mathutil.hpp"
#include <vector>

class Vector
{
private:
  int n; // The number of elements
  std::vector<double> v;
  void cleanUp(); // Deallocates memory
public:
  // Constructors and destructor
  Vector() : n(0) {} // Default constructor, zero length vector
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

inline Vector operator*(const Vector& v, const Matrix& mat)
{
  return lmultiply(v, mat);
}

inline Vector operator*(const Matrix& mat, const Vector& v)
{
  return rmultiply(mat, v);
}

#endif
