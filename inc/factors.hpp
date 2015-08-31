/*
 *    Purpose: Declare the signatures of several matrix factorisation routines
 *             and any associated procedures.
 * 
 *    DATE                AUTHOR               CHANGES
 *    =================================================================
 *    19/08/15            Robert Shaw          Original code.
 *    20/08/15            Robert Shaw          Added Householder.
 *    21/08/15            Robert Shaw          LU and Cholesky.
 *    22/08/15            Robert Shaw          Hessenberg added.
 */

#ifndef FACTORSHEADERDEF
#define FACTORSHEADERDEF

// Declare forward dependencies
class Matrix;
class Vector;
class Error;

// Declare the modified Gram-Schmidt procedure
// which takes a set of vectors in a full-rank
// matrix x, returning the orthogonalised vectors
// in a matrix q, and the transformation matrix r.
// Returns true if successful.
bool dgegs(const Matrix& x, Matrix& q, Matrix& r, const double& PRECISION); 

// Declare the Householder procedure
// giving the R matrix (y) of the QR factorisation of the matrix x
// and a set of reflection vectors, v, from which the q matrix 
// could be constructed. Returns true if successful 
bool dgehh(const Matrix& x, Matrix& y, Matrix& v);

// Implicity form product Qx using v from HH decomp
void implicitqx(const Matrix& v, Vector& x);

// Or Q(T)b instead
void implicitqtb(const Matrix& v, Vector& b);

// Return the full Q matrix from the HH decomp
Matrix explicitq(const Matrix& v);

// Get the LU decomposition of A by Gaussian Elimination with 
// partial pivoting. This actually computes PA = LU, putting L, U
// into the matrix B, and returning a vector of the row interchanges.
Vector dgelu(const Matrix& A, Matrix& B);

// Explicitly form the matrix P from the output of dgelu
Matrix explicitp(const Vector& p);

// Implicitly form Pb, where b is a vector, and P is the permutation
// matrix from the Gaussian eliminiation
void implicitpb(const Vector& p, Vector& b); 

// Compute the Cholesky factorisation A = R(T)R for a symmetric
// positive definite matrix, where R is an upper triangular matrix.
Matrix cholesky(const Matrix& A);

// Reduce a square matrix x into y in Hessenberg form, using Householder 
// reflections. It stores the reflectors in matrix v, which can be
// used implicitly later on. Returns true if successful.
bool hessenberg(const Matrix& x, Matrix& y, Matrix& v); 

// Procedures for computing and applying givens rotations:
// givens(a, b) will take scalars a, b and compute c = cos(t)
// and s=sin(t), returning them in the 2-vector [c, s].
// givens(A, G) will compute AG, while givens(G, A) will compute G(T)A
// lgivens(G, A) will compute GA.
// Here, G is the givens matrix, specified by [c, s], and A is a matrix
Vector givens(double a, double b, double PRECISION = 1e-12);
Matrix givens(const Vector& G, const Matrix& A, int i, int k);
Matrix givens(const Matrix& A, const Vector& G, int i, int k);
Matrix lgivens(const Vector& G, const Matrix& A, int i, int k);

// Explicitly form the givens matrix from the 2-vector g from givens(a, b)
// given matrix positions i and k, and dimension dim
Matrix explicitg(const Vector& g, int i, int k, int dim);

#endif
