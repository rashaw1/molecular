// Implements factors.hpp

#include "mvector.hpp"
#include "matrix.hpp"
#include "error.hpp"
#include "factors.hpp"
#include <iostream>
#include <cmath>

// Gram-Schmidt orthonormalises the columns of x, into the column 
// vectors of q, whilst generating the transformation matrix r,
// returning true if successful. Gives QR factorisation.
bool dgegs(const Matrix& x, Matrix& q, Matrix& r, const double& PRECISION)
{
  bool rVal = true; // Will only be changed to false if a problem occurs
  int dim = x.nrows(); // Assume that x is full rank
  // Make sure q and r are the appropriate sizes
  q.resize(dim, dim);
  r.resize(dim, dim);
  
  // Start procedure
  // Get columns of x as an array of vectors
  Vector* xa = new Vector[dim];
  for (int i = 0; i < dim; i++) {
    xa[i] = x.colAsVector(i);
  }
  // Define a placeholder array of vectors for q
  Vector* qa = new Vector[dim];
  r(0, 0) = pnorm(xa[0], 2); // Calculate 2-norm of x_0
  if (r(0, 0) < PRECISION) { // Singular norm - algorithm fails
    rVal = false;
  } else {
    qa[0] = (1.0/r(0,0))*xa[0]; // Normalise x_0
    // Begin main loop
    Vector y(dim); // Temporary placeholder
    for(int j = 1; j < dim; j++){
      y = xa[j]; 
      for(int i = 0; i < j; i++){
        r(i, j) = inner(qa[i], y); // qT*y
        y = y - r(i, j)*qa[i];
      }
      r(j, j) = pnorm(y, 2); // 2-norm of y
      if (r(j, j) < PRECISION) { 
        rVal = false; 
        break; 
      } else {
        // Normalise y into q
        qa[j] = (1.0/r(j, j))*y;
      }
    }
  }
  // Copy the vectors qa into q, and clean up memory
  for (int i = 0; i < dim; i++) {
    q.setCol(i, qa[i]);
  }
  delete[] xa;
  delete[] qa;
  return rVal;
}

// Householder factorisation algorithm, taking matrix x
// and transforming to the R matrix of QR factorisation
// returning a set of reflection vectors, v. Works on any 
// m x n matrix, with m >= n. 
bool dgehh(const Matrix& x, Matrix& y, Matrix& v)
{
  bool rVal = true; // Return value
  // Get dimensions of x
  int m = x.nrows();
  int n = x.ncols();
  // Make sure v and y are the same size as x
  v.resize(m, n); 
  y.resize(m, n);
  y = x;

  // Start main algorithm
  double value = 0.0; // Placeholder for norms
  for (int k = 0; k < n; k++){
    // Make temporary vector and matrix for subvector/matrix of x
    Vector column(m-k);
    Matrix subx(m-k, n-k);
    // Copy in values
    for (int i = k; i < m; i++){
      column[i-k] = y(i, k);
      for (int j = k; j < n; j++){
      	subx(i-k, j-k) = y(i, j);
      }
    }
    // Compute the reflector
    value = pnorm(column, 2);
    value = (column(0) < 0 ? -value : value);
    column[0] = column(0)+ value;
    value = pnorm(column, 2);
    column = (1.0/value)*column;

    // Transform the submatrix
    Vector temporary(n-k);
    temporary = column*subx;
    Matrix transx(m-k, n-k);
    transx = outer(column, temporary);
    subx = subx - 2.0*transx;

    // Transfer values to output matrices
    for (int i = 0; i < k; i++){
      v(i, k) = 0.0;
    }
    for (int i = k; i < m; i++){
      v(i, k) = column(i-k);
      for (int j = k; j < n; j++){
	y(i, j) = subx(i-k, j-k);
      }
    } 
  }
  return rVal;
}

// Take the matrix v from the HH decomposition and implicitly
// calculate the product of Q with a vector x, returned in x.
void implicitqx(const Matrix& v, Vector& x)
{
  int m = v.nrows();
  int n = v.ncols();
  // Begin main loop
  for (int k = n-1; k > -1; k--) {
    Vector temp(m-k); // Temporary vector of correct size
    Vector vk(m-k); 
    // Transfer values
    for (int i = 0; i < m-k; i++){
      temp[i] = x(i+k);
      vk[i] = v(i+k, k);
    }
    // Compute the product
    double dval = inner(vk, temp);
    vk = dval*vk;
    temp = temp - 2.0*vk;
    // Copy values back to output
    for (int i = 0; i < m-k; i++){
      x[i+k] = temp(i);
    }
  }
}

// Calculate the product Q(T)b of the Q matrix using the HH
// decomposition, with a vector b - implicitly, i.e. without
// ever forming Q
void implicitqtb(const Matrix& v, Vector& b)
{
  int m = v.nrows();
  int n = v.ncols();
  // Begin main loop
  for (int k = 0; k < n; k++) {
    Vector temp(m-k); // Temporary vector for subvector of b
    Vector vk(m-k); 
    // Copy values in
    for (int i = k; i < m; i++) {
      temp[i-k] = b(i);
      vk[i-k] = v(i, k);
    }
    double dval = inner(vk, temp);
    vk = dval*vk;
    temp = temp - 2.0*vk;
    // Copy values back to output
    for (int i = k; i < m; i++) {
      b[i] = temp(i-k);
    }
  }
}   

// Form explicitly the Q matrix from the HH decomposition
Matrix explicitq(const Matrix& v)
{
  int m = v.nrows();
  int n = v.ncols();
  Matrix rmat(m, n); // Return matrix
  // Do implicitqx for x = the identity vector in each dimension
  for (int i = 0; i < n; i++) {
    Vector temp(m, 0.0); // Temporary vector of all zeroes
    temp[i] = 1.0; // Turn into the identity
    implicitqx(v, temp); // Get the ith column of q
    rmat.setCol(i, temp); // Set the ith column of q 
  }
  return rmat;
}
  
// Gaussian elimination with partial pivoting LU decomposition
// assumes matrix is square (as nobody uses it otherwise!)
// Both L and U are stored in B.
// Returns a vector p with the order in which rows were interchanged
Vector dgelu(const Matrix& A, Matrix& B) 
{
  int dim = A.nrows(); // Assume square
  B = A; // Copy A into B
  Vector p(dim-1); // For returning the permutations at each step
  Matrix L(dim, dim, 0.0); // Matrix of all zeroes
  // Set L to be the identity
  for (int i = 0; i < dim; i++){
    L(i, i) = 1.0;
  }

  // Begin main algorithm
  for (int k = 0; k < dim-1; k++){
    // Find the pivot in column k
    int pivot = k;
    double testval = B(k, k);
    // Choose the pivot as the biggest (by absolute value)
    // element in column k
    for (int i = k+1; i < dim; i++){
      if(fabs(B(i, k)) > testval){
	pivot = i;
      }
    }
    // Interchange rows for submatrices of B and L
    B.swapRows(k, pivot, k); // Only swap latter part of row
    L.swapRows(k, pivot, 0, k); // Only swap front part of row
    // Store that at step k, row pivot was swapped
    p[k] = pivot;
    // Do the elimination step
    for (int j = k+1; j < dim; j++){
      L(j, k) = B(j, k)/B(k, k);
      for (int a = k; a < dim; a++){
	B(j, a) = B(j, a) - L(j, k)*B(k, a);
      }
    }
  }
  // B and L are now upper and lower triangular, respectively
  // and L has ones on the diagonal, so we store L in B
  for (int i = 1; i < dim; i++){ // Loop through rows
    for (int j = 0; j < i; j++) { // And portions of columns
      B(i, j) = L(i, j);
    }
  }
  return p;
}

// Explicitly form the permutation matrix from the output of dgelu
Matrix explicitp(const Vector& p)
{
  int dim = p.size()+1; // p will always m-1 values for an m x m matrix
  Matrix P(dim, dim, 0.0); // To return the permutation matrix in
  // Set P as the identity to start with
  for (int i = 0; i < dim; i++){
    P(i, i) = 1.0;
  }
  // Start construction
  for(int i = 0; i < dim-1; i++){
    Matrix temp(dim, dim, 0.0);
    // Need to maintain rows that weren't swapped
    for (int j = 0; j < dim; j++){
      if(j!= i && j!=p(i)){
	temp(j, j) = 1.0;
      }
    }
    temp(p(i), i) = 1.0; // rows i and p(i) were swapped at step i
    temp(i, p(i)) = 1.0; 
    P = temp*P;
  }
  return P;
}

// Implicity calculate Pb without ever forming P, 
// where P is the permutation matrix from the dgelu procedure
void implicitpb(const Vector& p, Vector& b)
{
  int dim = p.size();
  double temp = 0.0;
  for (int i = 0; i < dim; i++){
    // Swap the rows of b as specified by p
    temp = b(i);
    b[i] = b(p(i));
    b[p(i)] = temp;
  }
}

// Compute the Cholesky factorisation of a real symmetric
// positive definite matrix. Returns the upper triangular
// matrix R.
Matrix cholesky(const Matrix& A)
{
  int dim = A.nrows(); // Assume square
  Matrix R;
  R = A; // Initialise R
  // Reduce the elements of R symmetrically
  for (int k = 0; k < dim; k++){
    for (int j = k+1; j < dim; j++){
      for (int i = j; i < dim; i++){
	R(j, i) = R(j, i) - (R(k, j)/R(k, k))*R(k, i);
      }
    }
    double rootval = sqrt(R(k, k));
    for (int i = k; i < dim; i++){
      R(k, i) = R(k, i)/rootval;
    }
  }
  // Set sub-diagonal elements to be zero
  for (int i = 1; i < dim; i++){
    for (int j = 0; j < i; j++){
      R(i, j) = 0.0;
    }
  }
  return R;
}

// Decompose the square matrix x into hessenberg form in y,
// giving the householder reflectors in v. Returns true if successful.
bool hessenberg(const Matrix& x, Matrix& y, Matrix& v)
{
  bool rval = true;
  int dim = x.nrows();
  if (x.isTriangular()){ // No need to reduce it
    y = x;
    v.assign(dim, dim, 0.0);
    for (int i = 0; i < dim; i++) { v(i, i) = 1.0; }
  }
  else if (x.isSquare()){ // Must be square to proceed    
    // Make sure v, y are the right size
    y = x; // Initialise to input matrix
    v.assign(dim, dim, 0.0); // Make all zeroes 
    // Begin main loop
    for (int k = 0; k < dim-2; k++){
      bool reduce = false;
      int check = k+2;
      while(check < dim && !reduce){ reduce = (fabs(y(check, k)) > 1e-14 ? true : false); check++; }
      if(reduce){
      // Declare temporary vectors
      Vector xk(dim-k-1); 
      Vector vk(dim-k-1);
      // Copy in the values needed
      for (int i = k+1; i < dim; i++){
		xk[i-k-1] = y(i, k);
      }
      // Compute the reflector
      // Choose the sign that maxmises distance of reflected vector
      double val = (xk(0) < 0 ? -1.0*pnorm(xk, 2) : pnorm(xk, 2));
      vk = xk;
      vk[0] += val; // Correct leading value
      vk  = (1.0/pnorm(vk, 2))*vk; // Normalise
      // Compute the changes to y
      Matrix temp1(dim-k-1, dim-k);
      // Copy values needed into temp matrix
      for (int i = k+1; i < dim; i++){
	for (int j = k; j < dim; j++){
	  temp1(i-k-1, j-k) = y(i, j);
	}
      }
      Vector temp2(dim-k-1);
      temp2 = vk*temp1;
      temp1 = temp1 - 2.0*outer(vk, temp2);
      // Copy values back in to y
      for (int i = k+1; i < dim; i++){
	for (int j = k; j < dim; j++){
	  y(i, j) = temp1(i-k-1, j-k);
	}
      }
      // Repeat on other side
      temp1.resize(dim, dim-k-1);
      // Copy values into temp matrix
      for (int i = 0; i < dim; i++){
	for (int j = k+1; j < dim; j++){
	  temp1(i, j-k-1) = y(i, j);
	}
      }
      temp2 = temp1*vk;
      temp1 = temp1 - 2.0*outer(temp2, vk);
      // Copy values back into y
      for (int i = 0; i < dim; i++){
	for (int j = k+1; j < dim; j++){
	  y(i, j) = temp1(i, j-k-1);
	}
      }
      // Copy reflection vector into v
      for (int i = k+1; i < dim; i++){
	v(i, k+1) = vk(i-k-1);
      }
    }
    }
  } else {
    rval = false; // Algorithm failed
  }
  return rval;
}

// Compute and apply givens rotations
Vector givens(double a, double b, double PRECISION)
{
  Vector g(2); // Return vector, [c, s]
  double c, s;
  if (fabs(b) < PRECISION){
    c = 1.0; s=0.0;
  } else {
    double tau;
    if (fabs(b) - fabs(a) > PRECISION){
      tau = -1.0*a/b;
      s = 1.0/sqrt(1+tau*tau);
      c = s*tau;
    } else {
      tau = -1.0*b/a;
      c = 1.0/sqrt(1+tau*tau);
      s = c*tau;
    }
  }
  g[0] = c; g[1] = s;
  return g;
}

// G is the givens rotation G(i, k, t), represented by the vector [c, s] and
// the positions i, k.
Matrix givens(const Vector& G, const Matrix& A, int i, int k)
{
  Matrix GA; // Will return the product G(T)A
  GA = A;
  double tau1, tau2;
  for (int j = 0; j < A.ncols(); j++){
    tau1 = A(i, j);
    tau2 = A(k, j);
    // Only two rows are affected
    GA(i, j) = G(0)*tau1 - G(1)*tau2;
    GA(k, j) = G(1)*tau1 + G(0)*tau2;
  }
  return GA;
}

Matrix lgivens(const Vector& G, const Matrix& A, int i, int k){
  Matrix GA; // Will return the product GA
  GA = A;
  double tau1, tau2;
  for (int j = 0; j < A.ncols(); j++){
    tau1 = A(i, j);
    tau2 = A(k, j);
    // Only two rows affected
    GA(i, j) = G(0)*tau1 - G(1)*tau2;
    GA(k, j) = G(0)*tau2 + G(1)*tau1;
  }
  return GA;
}

Matrix givens(const Matrix& A, const Vector& G, int i, int k){
  Matrix AG; // Will return the product AG
  AG = A;
  double tau1, tau2;
  for (int j = 0; j < A.nrows(); j++){
    tau1 = A(j, i);
    tau2 = A(j, k);
    // Only two columns affected
    AG(j, i) = G(0)*tau1 - G(1)*tau2;
    AG(j, k) = G(1)*tau1 + G(0)*tau2;
  }
  return AG;
}

// Explicitly form a givens matrix
Matrix explicitg(const Vector& g, int i, int k, int dim)
{
  Matrix G(dim, dim, 0.0); // Make a dim x dim matrix of zeroes
  // Turn into the identity matrix, except with c at positions
  // (i, i) and (k, k)
  for (int j = 0; j < dim; j++){
    if(j == i || j == k){
      G(j, j) = g(0);
    } else {
      G(j, j) = 1.0;
    }
  }
  // Set the s elements
  G(i, k) = g(1);
  G(k, i) = -1.0*g(1);
  return G;
}
  
