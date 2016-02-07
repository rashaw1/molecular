/*
 *     PURPOSE: defines and implements class Matrix, a matrix with members
 *              of any class T that has defined upon it the usual
 *              arithemetic operations.
 *
 *     NOTE: implemented in header file due to templating - avoids some
 *           potential linking issues.
 *
 *     DATE             AUTHOR                CHANGES
 *     =====================================================================
 *     13/08/15         Robert Shaw           Original code
 *     14/08/15         Robert Shaw           Added error throwing
 *     20/08/15         Robert Shaw           Changed approach to matrix-
 *                                            matrix multiplication.
 */

#ifndef MATRIXHEADERDEF
#define MATRIXHEADERDEF

// Declare forward dependencies

#include "error.hpp"
#include "mvector.hpp"
#include <vector>

class Matrix
{
private:
  int rows, cols; // No. of rows and columns of the matrix
  std::vector<Vector> arr;
  void cleanUp(); // Utility function for memory deallocation
public:
  // Constructors and destructor
	Matrix(); // Default constructor
  Matrix(int m, int n); // Declare an m x n matrix
  Matrix(int m, int n, const double& a); // Declare m x n matrix, all entries = a
  Matrix(int m, int n, const Vector& a); // Matrix of m row copies of n-vector a
  Matrix(const Matrix& other); // Copy constructor
  ~Matrix(); // Destructor
  // Accessors
  int nrows() const { return rows; } // Returns no. of rows
  int ncols() const { return cols; } // Returns no. of cols  
  Vector rowAsVector(int r) const; //Returns row r as a vector
  Vector colAsVector(int c) const; //Returns col c as a vector
  void setRow(int r, const Vector& u); // Sets row r to be the vector u
  void setCol(int c, const Vector& u); // Sets col c to be the vector u
  // Shaping functions
  void resize(int m, int n); // Resize to empty m x n matrix
  void assign(int m, int n, const double& a); // Resize, setting all entries to a
  void removeRow(int r); // Remove row r
  void removeCol(int c); // Remove column c
  void swapRows(int i, int j, int start, int end); // Swap (sections of )rows i and j
  void swapCols(int i, int j, int start, int end); // Swap (sections of )cols i and j
  void swapRows(int i, int j, int start = 0); // Assumed start/end points
  void swapCols(int i, int j, int start = 0); // Assumed start/end points
  // Overloaded operators
  double& operator[](int i); // Return pointer to first element of row i
  double& operator()(int i, int j); // Return pointer to element ij
  double operator()(int i, int j) const; // Return by value element ij
  Matrix& operator=(const Matrix& other); 
  // Unary operators
  Matrix operator+() const; 
  Matrix operator-() const;
  // Binary operators
  Matrix operator+(const Matrix& other) const;
  Matrix operator-(const Matrix& other) const;
  Matrix& operator*=(const double& scalar) { return *this; } // Scalar multiplication
  Matrix operator*(const Matrix& other) const; // Matrix x matrix
  // Intrinsic functions
  Matrix transpose() const; // Return the transpose of the matrix
  double trace() const;
  void print(double PRECISION = 1e-12) const; // Pretty prints the matrix to primary ostream 
  bool isSymmetric() const; // Determines whether the matrix is symmetric
  bool isSquare() const { return ( rows == cols ); }  // Determines whether the matrix is square
  bool isTriangular(bool upper = true) const; // Determines whether it is upper/lower triangular
  // Friend functions
  friend double pnorm(const Matrix& m, int p); // The induced matrix p-norm
  friend double fnorm(const Matrix& m); // Calculate the Frobenius norm
};

// Scalar multiplication

inline Matrix operator*(const double& scalar, const Matrix& m) 
{
  // Make return matrix of correct size
  Matrix rMat(m.nrows(), m.ncols());
  // Multiply elements by scalar
  for (int i = 0; i < m.nrows(); i++){
    for (int j = 0; j < m.ncols(); j++){
      rMat(i, j) = m(i, j)*scalar;
    }
  }
  return rMat;
}

inline Matrix operator*(const Matrix& m, const double& scalar) 
{
  // Make return matrix of correct size
  Matrix rMat(m.nrows(), m.ncols());
  // Multiply elements by scalar
  for (int i = 0; i < m.nrows(); i++){
    for (int j = 0; j < m.ncols(); j++){
      rMat(i, j) = m(i, j)*scalar;
    }
  }
  return rMat;
}

#endif
