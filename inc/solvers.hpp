/*
 *   Purpose: To define a number of solvers, particularly for linear
 *            equations of the form Ax = b.
 * 
 *   DATE             AUTHOR            CHANGES
 *   ===============================================================
 *   21/08/15         Robert Shaw       Original code
 *   22/08/15         Robert Shaw       Iterative eigenv's  added.
 *   23/08/15         Robert Shaw       Symqr with implicit shifts.
 */

#ifndef SOLVERSHEADERDEF
#define SOLVERSHEADERDEF

// Declare forward dependencies
class Vector;
class Matrix;
class Error;

// The basic back-substitution routine, which is used in pretty much
// every other solver. R is an upper triangular matrix, y is the 
// rhs vector - i.e., solving Rx = y, where x is the solution returned
Vector backsub(const Matrix& R, const Vector& y);

// Solve Ax = b by QR factorisation (using Householder algorithm)
// construct y = Q(T)b, by implicitqtb, then solve Rx = y by 
// back substitution. 
// First instance performs decomposition, second takes v from 
// already performed decomposition.
Vector qrsolve(const Matrix& A, const Vector& b); 
Vector qrsolve(const Matrix& R, const Matrix& v, const Vector& b);

// Solve the full-rank least squares problem by QR factorisation
// Ax = y with A being an m x n matrix, m > n
Vector qrsquares(const Matrix& A, const Vector& b);

// Solve the square Ax = b problem by LU decomposition, i.e 
// Gaussian elimination with partial pivoting. 
// First instance does decomposition, second instance
// takes already formed decomposition.
Vector lusolve(const Matrix& A, const Vector& b);
Vector lusolve(const Matrix& B, const Vector& p, const Vector& b);

// Solve the square, symmetric positive definite sytem Ax = b 
// using cholesky factorisation
// Second instance uses already formed factorisation - note
// that the arguments are reversed, to distinguish the two
Vector choleskysolve(const Matrix& A, const Vector& b);
Vector choleskysolve(const Vector& b, const Matrix& R);


// Use the Penrose-Moore pseudo-inverse to solve a singular system
Vector stableSolve(const Matrix& A, const Vector& b);

// Use the power iteration algorithm to find an eigenvalue -
// specifically the eigenvalue with largest absolute value -
// of A given a (normalised) vector, v. The eigenvector is
// left in v. Returns the eigenvalue.
// Cuts off when PRECISION is reached, or MAXITER iterations done
double poweriter(const Matrix& A, Vector& v, double PRECISION = 1e-12, int MAXITER = 100);

// Use the inverse power iteration method to find the eigenvalue
// (and corresponding eigenvector) of A nearest to u, given a vector v.
// Returns the eigenvalue.
double inverseiter(const Matrix& A, Vector& v, double u, double PRECISION = 1e-12, int MAXITER=100);

// Use the Rayleigh Iteration Algorithm to get approximations for an eigenvec/val pair
// of a matrix A. Returns eigenvalue, stores vector in v.
double rayleigh(const Matrix& A, Vector& v, double l0, double PRECISION = 1e-8, int MAXITER = 50);

// Use the QR algorithm with shifts to find the approximate eigenvalues
// of a matrix A. The values are returned in the vector vals.
// The second instance is for when vectors are wanted as well.
// Will return true if successful.
bool qrshift(const Matrix& A, Vector& vals, double PRECISION = 1e-12, int MAXITER = 100);
bool qrshift(const Matrix& A, Vector& vals, Matrix& vecs, double PRECISION = 1e-12, int MAXITER=100);

// Implcitshift step needed for symqr                                           
Matrix implicitshift(Matrix& T, double PRECISION);

// The QR algorithm for a real-symmetric matrix using implicit shifts is more efficient
// than the above alternative
bool symqr(const Matrix& A, Vector& vals, double PRECISION = 1e-12);
bool symqr(const Matrix& A, Vector& vals, Matrix& vecs, double PRECISION, bool HESS = false);


// Utility functions for symeig that pack and unpack matrices              
void splitmatrix(const Matrix& B, Matrix& b1, Matrix& b2, int i);
void joinmatrix(Vector& vals, const Vector& vals1, const Vector& vals2,
		Matrix& vecs, const Matrix& vecs1, const Matrix& vecs2, int i);

/* UNDER CONSTRUCTION
// Use a divide and conquer algorithm to find the eigenvalues and vectors
// of a real symmetric matrix A. This is generally the fastest method
// when a complete set of both values and vectors are needed.
// Will return values in vector vals, vectors in the columns of matrix
// vecs, and returns true if successful.
bool symeig(const Matrix& A, Vector& vals, Matrix& vecs, double PRECISION = 1e-12, bool tridiag = false, int MINDAC = 10);

// Find zero of secular equation for symeig
double findzero(Vector& d, const Vector& v, int i, double PRECISION = 1e-12);

// Do the rank1 update part of the symeig procedure
void diagupdate(Vector& D, double bm, Vector& z, Vector& vals, Matrix& Q, double PRECISION);
*/

#endif
