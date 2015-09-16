// Implements solvers.hpp

#include "factors.hpp"
#include "solvers.hpp"
#include "mvector.hpp"
#include "matrix.hpp"
#include "error.hpp"
#include <cmath>
#include <iostream>

// Back substitution of the triangular system Rx = y
Vector backsub(const Matrix& R, const Vector& y)
{
  int dim = R.nrows(); // Triangular matrices must be square
  Vector rvec(dim); // Return vector
  // Begin algorithm
  double sum; // For storing the summation step in each loop
  for (int k = dim-1; k > -1; k--){
    // Add contributions from all previously determined entries
    sum = 0.0;
    for (int i = k+1; i < dim; i++){
      sum += rvec(i)*R(k, i);
    }
    // Determine final value
    rvec[k] = (y(k) - sum)/R(k, k);
  }
  return rvec;
}

// Solve the linear system Ax = b, using Householder-based QR
// factorisation, and backsub. Returns the solution vector x.
// Assumes A is square and nonsingular.
Vector qrsolve(const Matrix& A, const Vector& b)
{
  int dim = A.nrows();
  Vector x(dim); // For returning the answer
  Matrix v; Matrix r; // For the Householder algorithm
  // dgehh resizes v and r for us, so they
  // need only be declared here as null vectors.

  // QR Factorise A into v and r
  if(dgehh(A, r, v)){
    // Construct y = Q(T)b (x is y to save space)
    x = b; // Copy b into x
    implicitqtb(v, x);
    // Solve Rx = y by backsubstitution
    x  = backsub(r, x);
  } else {
    throw( Error("QRFACT", "Householder triangularisation unsuccesful.") );
  }
  return x;
}

// Second instance does the same as above, but is given the matrices R and v as
// arguments, so repeated decompositions can be avoided.
Vector qrsolve(const Matrix& R, const Matrix& v, const Vector& b)
{
  Vector x(b.size()); // Solution vector
  x = b; // Initialise x to b
  implicitqtb(v, x); // Calculate Q(T)b implicitly
  // Solve Rx = y = Q(T)b by backsub
  x = backsub(R, x);
  return x;
}

// Use QR factorisation to solve the full-rank least-squares problem
Vector qrsquares(const Matrix& A, const Vector& b)
{
  int m = A.nrows();
  int n = A.ncols();
  Vector x(n); // Solution vector
  if ( m < n ) {
    throw( Error("QRSQRS", "Least squares problem is rank-deficient.") );
  } else {
    Matrix r; Matrix v;
    // Get the QR factorisation
    if(dgehh(A, r, v)){
      x = b; // Store b in x for convenience
      implicitqtb(v, x); // Calculate Q(T)b
      // Reduce the R matrix and Q(T)b vector
      for (int i = m-1; i > n-1; i--){
	r.removeRow(i);
      }
      x.resizeCopy(n);
      x = backsub(r, x);
    }
  }
  return x;
}

// Use the LU decomposition (Gaussian elimination with partial pivoting)
// to solve the system Ax = b. First, PA = LU is formed, and so we must 
// solve PAx = Pb. Pb is formed by implicitpb, then we can solve LUx = Pb
// first by forward substitution Ly = Pb for  y = Ux, then solve Ux = y
// by back substitution for x
Vector lusolve(const Matrix& A, const Vector& b)
{
  int dim = A.nrows(); // Assume square
  Vector x(dim); // Will contain the solution
  // First, LU decompose A
  Vector p(dim-1);
  Matrix B(dim, dim);
  p = dgelu(A, B);
  // Calculate Pb implicitly
  x = b;
  implicitpb(p, x);
  // Solve Ly = Pb by forward substitution
  // remembering diagonal of L is all ones
  for (int i = 1; i < dim; i++){
    double sum = 0.0;
    for (int j = 0; j < i; j++){
      sum += B(i, j)*x[j];
    }
    x[i] = x[i] - sum;
  }
  // Now we need to solve Ux = y for x
  // First remove L from B, to get U
  for(int i = 1; i<dim; i++){
    for (int j = 0; j < i; j++){
      B(i, j) = 0.0;
    }
  }
  // Backsubstitute to get x
  x = backsub(B, x);
  return x;
}

// Do the same as above, but with already formed LU decomp so as to
// avoid the need for repeated decompositions
Vector lusolve(const Matrix& B, const Vector& p, const Vector& b)
{
  int dim = b.size();
  Vector x(dim); // Solution vector
  x = b;
  implicitpb(p, x); // Calculate Pb implicitly
  // Solve Ly = Pb by forward substitution
  //remembering diagonal of L is all ones
  for (int i = 1; i < dim; i++){
    double sum = 0.0;
    for (int j = 0; j < i; j++){
      sum += B(i, j)*x[j];
    }
    x[i] = x[i] - sum;
  }
  // Now we need to solve Ux = y for x                     
  // First remove L from B, to get U
  Matrix U;
  U = B;
  for(int i = 1; i<dim; i++){
    for (int j = 0; j < i; j++){
      U(i, j) = 0.0;
    }
  }
  // Then solve by backsub
  x = backsub(U, x);
  return x;
}


// Solve the linear system Ax = b, where A is real symmetric
// positive definite, using Cholesky decomposition
// The algorithm is A = R(T)R by decomposition, so we solve
// R(T)Rx = b - so solve the lower triangular R(T)y = b
// by forward substitution, then solve the upper triangular
// Rx = y by back substitution.
Vector choleskysolve(const Matrix& A, const Vector& b)
{
  int dim = A.nrows(); // Assume square
  Vector x(dim); // Solution vector
  Matrix R;
  R = cholesky(A); // Get the upper triangular matrix R

  // Do the forward substitution for y
  x = b; // Initialise x
  for (int i = 0; i < dim; i++){
    double sum = 0.0;
    for (int j = 0; j < i; j++){
      sum += R(j, i)*x[j]; // We need the transpose of R
    }
    x[i] = x[i] - sum;
    x[i] = x[i]/R(i, i); // Normalise
  }
  // Now do the back substitution
  x = backsub(R, x);
  return x;
}

// Do the same as above where the decomposition is already given -
// note the arguments are reversed
Vector choleskysolve(const Vector& b, const Matrix& R)
{
  int dim = b.size();
  Vector x(dim); // Solution vector
  x = b;
  // Do the forward substitution for y 
  for (int i = 0; i < dim; i++){
    double sum = 0.0;
    for (int j = 0; j < i; j++){
      sum += R(j, i)*x[j]; // We need the transpose of R           
    }
    x[i] = x[i] - sum;
    x[i] = x[i]/R(i, i); // Normalise                                                   
  }
  // Now do the back substitution                                         
  x = backsub(R, x);
  return x;
}

Vector stableSolve(const Matrix& A, const Vector& b)
{
  Vector x(b.size());
  return x;
}

// The power iteration method for finding the largest absolute 
// valued eigenvalue of a matrix A, given a normalised vector v.
double poweriter(const Matrix& A, Vector& v, double PRECISION, int MAXITER)
{
  Vector w; // Temporary vector for intermediate steps
  int dim = v.size();
  double norm = pnorm(v, 2);
  double lambda = 0.0; // The eigenvalue 
  if(norm - 1.0 > PRECISION) { // Normalise if not already
    v = (1.0/norm)*v;
  }
  // Begin loop
  double dist = 1.0; // Track distance between eigenvalue at each iter
  double err = 1.0; // Track error - max of dist and norm
  double oldlambda = 0.0; // Store the old lambda value for error calc
  int iter = 0; // Track number of iterations
  Vector oldv;
  while(err > PRECISION && iter < MAXITER){
    oldv = v; // Store previous vector
    w = A*v; // Calculat Av
    v = w;
    w.sort(); // Sort w without changing v
    // The best guess for the eigenvalue is the largest by absolute value
    // member of w
    lambda = (fabs(w(0)) > fabs(w(dim-1)) ? w(0) : w(dim-1));
    v = (1.0/lambda)*v;
    dist = fabs(lambda - oldlambda);
    norm = pnorm(v-oldv, 2);
    err = (norm < dist ? dist : norm);
    oldlambda = lambda;
    iter++;
  }
  return lambda;
}

// Do the inverse power iteration to find the eigenvalue closest to u, and the
// corresponding eigenvector (stored in v). Uses lusolve.
double inverseiter(const Matrix& A, Vector& v, double u, double PRECISION, int MAXITER)
{
  Vector w; // Temporary vector for intermediate steps                              
  int dim = v.size();
  double norm = pnorm(v, 2);
  double lambda = 0.0; // The eigenvalue                                   
  if(norm - 1.0 > PRECISION) { // Normalise if not already                     
    v = (1.0/norm)*v;
  }
  // Form the matrix A-uI and LU decompose
  Matrix X(dim, dim);
  X = A;
  for (int i = 0; i < dim; i++){
    X(i, i) -= u;
  }
  Matrix B; Vector p; // Will store the LU/cholesky decomposition
  p = dgelu(X, B);

  // Begin loop                                                                   
  double dist = 1.0; // Track distance between eigenvalue at each iter                  
  double err = 1.0; // Track error - max of dist and norm                            
  double oldlambda = 0.0; // Store the old lambda value for error calc              
  int iter = 0; // Track number of iterations                                 
  Vector oldv;
  while(err > PRECISION && iter < MAXITER){
    oldv = v; // Store previous vector                           
    // Solve the system of equations
    w = lusolve(B, p, v);
    v = w;
    w.sort();
    lambda = (fabs(w(0)) > fabs(w(dim-1)) ? w(0) : w(dim-1));
    v = (1.0/lambda)*v;
    lambda = (1.0/lambda) + u;
    dist = fabs(lambda - oldlambda);
    norm = pnorm(v-oldv, 2)/(pnorm(v, 2));
    err = (norm < dist ? dist : norm);
    oldlambda = lambda;
    iter++;
  }
  v = (1.0/pnorm(v, 2))*v;
  return lambda;
}  

// Rayleigh quotient iteration to get fast convergence to an eigenvalue
// eigenvector pair. Takes the matrix A for which eigenv's are needed,
// a starting guess eigenvec v, a starting guess value l0, and cutoff criteria.
// Returns eigenvalue, stores vector in v.
double rayleigh(const Matrix& A, Vector& v, double l0, double PRECISION, int MAXITER)
{
  double lambda = 0.0;
  double mu = l0;
  v = (1.0/pnorm(v, 2))*v;
  // Form A - mu*I
  Matrix X;
  X = A;
  for (int i = 0; i < A.nrows(); i++){
    X(i, i) = X(i, i) - mu;
  }
  // Solve Xy = v for y
  Vector y;
  y = lusolve(X, v);
  lambda = inner(y, v);
  mu = mu + (1.0/lambda);
  double err = pnorm(y-lambda*v, 2)/(pnorm(y, 2));
  int iter = 0;
  while (err > PRECISION && iter < MAXITER){
    v = (1.0/pnorm(y, 2))*y;
    // Form A - mu*I
    X = A;
    for (int i = 0; i < A.nrows(); i++){
      X(i, i) = X(i, i) - mu;
    }
    // Solve Xy = v
    y = lusolve(X, v);
    lambda = inner(y, v);
    mu = mu + (1.0/lambda);
    err = pnorm(y-lambda*v, 2)/pnorm(y, 2);
    iter++;
  }
  return mu;
}

// QR Algorithm with shifts to calculate the eigenvalues of a matrix A.
// Assumes real symmetric matrix.
bool qrshift(const Matrix& A, Vector& vals, double PRECISION, int MAXITER)
{
  bool rval = true;
  int n = A.nrows(); // Must be square
  vals.resize(n);
  Matrix v;
  Matrix vectemp;
  // Reduce A to tridiagonal form
  if(hessenberg(A, vectemp, v)){
    // Begin main loop
    for (int m = n-1; m > 0; m--){
      // Proceed until subdiagonal element is essentially zero
      int iter = 0;
      while(fabs(vectemp(m-1, m)) > PRECISION && iter < MAXITER){
	// Form vecs - vec(m, m)*I = temp
	Matrix temp1;
	temp1 = vectemp;
	for (int i = 0; i<m+1; i++){
	  temp1(i, i) = temp1(i, i) - temp1(m, m);
	}
	// Do householder qr decomp
	Matrix temp2;
	if (dgehh(temp1, temp2, v)){
	  // Form R*Q
	  temp1 = explicitq(v);
	  temp1 = temp2*temp1;
	  // Add in vectemp(m, m) * I
	  for (int i = 0; i < m+1; i++){
	    temp1(i, i) = temp1(i, i) + vectemp(m, m);
	  }
	  vectemp = temp1;
	} else {
	  rval = false; // Something went wrong
	}
	iter++;
      }
      // Record eigenvalue m
      vals[m] = vectemp(m, m);
      // Deflate vectemp
      vectemp.removeRow(m);
      vectemp.removeCol(m);
    }
    // Record eigenvalue 0
    vals[0] = vectemp(0, 0);
  } else {
    rval = false; // Something went wrong with hessenberg
  }
  return rval;
}

// Same as above, but for when vectors are wanted as well
bool qrshift(const Matrix& A, Vector& vals, Matrix& vecs, double PRECISION, int MAXITER)
{
  bool rval = true;
  if(qrshift(A, vals, PRECISION, MAXITER)){
    // Find eigenvectors using inverse iteration, as we have a good guess for eigenvalue
    double temp;
    int dim = vals.size();
    Vector v(dim);
    for (int i = 0; i<dim; i++){
      v[i] = 1.0;
    }
    vecs.resize(dim, dim);
    for(int i = 0; i < dim; i++){
      temp = inverseiter(A, v, vals(i), PRECISION, MAXITER);
      vecs.setCol(i, v);
    }
  } else {
    rval = false;
  }
  return rval;
}

// The real symmetric case is more efficiently solved by using implicit shifts
// as in the following implementation:
bool symqr(const Matrix& A, Vector& vals, double PRECISION)
{

  bool rval = true;
  int dim = A.nrows(); // It's square
  vals.resize(dim);
  Matrix B; Matrix q;
  // Tridiagonalise
  if(hessenberg(A, B, q)){
    int flag = 0;
    while (flag < dim-1){
      // Reduce B
      for (int i = 0; i < dim-1; i++){
	// Set elements to zero if tiny
	if (fabs(B(i, i+1)) < PRECISION){
	  B(i, i+1) = B(i+1, i) = 0.0;
	}
      }
      // See if diagonal matrix in top left
      int p = 0;
      while (p < dim - 1){
	if(B(p, p+1) == 0.0){
	  p++;
	} else {
	  break;
	}
      }
      // See if diagonal matrix in bottom right
      int q = dim-1; 
      while(q > p){
	if(B(q, q-1) == 0.0){
	  q--;
	} else {
	  break;
	}
      }
      // p is now the point at which the unreduced tridiagonal matrix begins
      // and q is the point at which the lower right diagonal matrix begins
      
      // We now perform the shift on the unreduced matrix
      if (p!=q && flag < dim-1) { // If p = q, then the matrix is diagonal already!
	Matrix D(q-p+1, q-p+1, 0.0);
	// Copy values in
	for (int i = p; i < q; i++){
	  D(i-p, i-p) = B(i, i);
	  D(i-p+1, i-p) = D(i-p, i-p+1) = B(i, i+1);
	}
	D(q-p, q-p) = B(q-p, q-p);
	// Do the implicit shift step, getting the transformation matrix Z
	Matrix Z;
	Z = implicitshift(D, PRECISION);
	// Recompute B
	for (int i = p; i < q; i++){
	  B(i, i) = D(i-p, i-p);
	  B(i, i+1) = B(i+1, i) = D(i-p+1, i-p);
	}
	B(q-p, q-p) = D(q-p, q-p);
      }
      flag = p;
    }
    // Copy eigenvalues from diagonal of B
    for (int i = 0; i < dim; i++){
      vals[i] = B(i, i);
    }
  } else { // Hessenberg failed
    rval = false;
  }
  return rval;
}

// Same as above, but computes eigenvectors as well
bool symqr(const Matrix& A, Vector& vals, Matrix& vecs, double PRECISION, bool HESS)
{
  
  bool rval = true;
  int dim = A.nrows(); // It's square
  vals.resize(dim);
  vecs.assign(dim, dim, 0.0);
  if (dim == 1) { vals[0] = A(0, 0); vecs(0, 0) = 1.0; }
  else if (dim == 2) { 
    // Calculate analytically
    double b = A(0, 0) + A(1, 1);
    double ac = A(0, 0)*A(1,1) - A(0, 1)*A(1, 0);
    double b24ac = std::sqrt(b*b - 4*ac);
    vals[0] = 0.5*(b - b24ac);
    vals[1] = 0.5*(b + b24ac);
    double mult = A(0, 1)/(vals(0) - A(0,0));
    vecs(1, 0) = 1/(std::sqrt(1+mult*mult));
    vecs(0, 0) = mult*vecs(1, 0);
    mult = A(0, 1)/(vals(1) - A(0, 0));
    vecs(1, 1) = 1/(std::sqrt(1+mult*mult));
    vecs(0, 1) = mult*vecs(1, 1);
  } else {
  Matrix B; Matrix q;
  // Tridiagonalise 
  if (HESS) {
    B = A; 
    A.print();
    B.print();
    for (int i = 0; i < dim; i++) vecs(i, i) = 1.0;
  } else {
    rval = hessenberg(A, B, q);
    vecs = explicitq(q);
  }
  if(rval){
    
    int flag = 0;
    while (flag < dim-1){
      // Reduce B
      for (int i = 0; i < dim-1; i++){
	// Set elements to zero if tiny
	if (fabs(B(i, i+1)) < PRECISION){
	  B(i, i+1) = B(i+1, i) = 0.0;
	}
      }
      // See if diagonal matrix in top left
      int p = 0;
      while (p < dim - 1){
	if(B(p, p+1) == 0.0){
	  p++;
	} else {
	  break;
	}
      }
      // See if diagonal matrix in bottom right
      int q = dim-1; 
      while(q > p){
	if(B(q, q-1) == 0.0){
	  q--;
	} else {
	  break;
	}
      }
      // p is now the point at which the unreduced tridiagonal matrix begins
      // and q is the point at which the lower right diagonal matrix begins
      
      // We now perform the shift on the unreduced matrix
      if (p!=q && flag < dim-1) { // If p = q, then the matrix is diagonal already!
	Matrix D(q-p+1, q-p+1, 0.0);
	// Copy values in
	for (int i = p; i < q; i++){
	  D(i-p, i-p) = B(i, i);
	  D(i-p+1, i-p) = D(i-p, i-p+1) = B(i, i+1);
	}
	D(q-p, q-p) = B(q-p, q-p);

	// Check for more zeroes on superdiagonal
	int flag2 = 0;
	int bpoint = 0;
	while (flag2 < q-p) {
	  if(D(flag2, flag2+1) == 0 ){
	    bpoint = flag2;
	    flag2 = q-p;
	  } else {
	    flag2++;
	  }
	}    
	 
	// If block diagonal, break it up and call recursively
	if (bpoint != 0){
	  D.print();
	  Matrix D1(bpoint+1, bpoint+1, 0.0);
	  Matrix D2(q-p-bpoint, q-p-bpoint, 0.0);
	  // Copy values in
	  for (int i = 0; i < bpoint; i++){
	    D1(i, i) = D(i, i);
	    D1(i, i+1) = D1(i+1, i) = D(i, i+1);
	  }
	  D1(bpoint, bpoint) = D(bpoint, bpoint);
	  for (int i = bpoint+1; i < q-p; i++){
	    D2(i-bpoint-1, i-bpoint-1) = D(i, i);
	    D2(i-bpoint-1, i-bpoint) = D2(i-bpoint, i-bpoint-1) = D(i, i+1);
	  }
	  D2(q-p-bpoint-1, q-p-bpoint-1) = D(q-p, q-p);
	  D1.print(); D2.print();
	  // Call symqr recursively
	  Vector vals1; Vector vals2;
	  Matrix D1P; Matrix D2P;
	  rval = symqr(D1, vals1, D1P, PRECISION, true) && rval;
	  rval = symqr(D2, vals2, D2P, PRECISION, true) && rval;
	  
	  // Construct the solution
	  for (int i = 0; i < bpoint+1; i++) {
	    B(i+p, i+p) = vals1(i);
	    for (int j = 0; j < bpoint+1; j++)
	      vecs(i+p, j+p) = D1P(i, j);
	  }
	  for (int i = bpoint+1; i < q-p+1; i++){
	    B(i+p, i+p) = vals2(i-bpoint-1);
	    for (int j = bpoint+1; j < q-p+1; j++)
	      vecs(i+p, j+p) = D2P(i-bpoint-1, j-bpoint-1);
	  }
	  flag = dim;
	} else {
	  // Do the implicit shift step, getting the transformation matrix Z
	  Matrix Z;
	  Z = implicitshift(D, PRECISION);
	  // Recompute B
	  for (int i = p; i < q; i++){
	    B(i, i) = D(i-p, i-p);
	    B(i, i+1) = B(i+1, i) = D(i-p+1, i-p);
	  }
	  B(q-p, q-p) = D(q-p, q-p);
	  // Recompute Q
	  D.assign(dim, dim, 0.0);
	  for (int i = 0; i < p; i++) { D(i, i) = 1.0; }
	  for (int i = q+1; i < dim; i++) { D(i, i) = 1.0;}
	  for (int i = p; i < q+1; i++){
	    for (int j = p; j < q+1; j++){
	      D(i, j) = Z(i-p, j-p);
	    }
	  }
	  vecs = vecs*D;
	  flag = p;
	} 
      } else {
	flag = p;
      }
      
    }
    // Copy eigenvalues from diagonal of B
    for (int i = 0; i < dim; i++){
      vals[i] = B(i, i);
    }
  } else { // Hessenberg failed
    rval = false;
  }
  }
  return rval;
}

// This does the implicit symmetric QR step with Wilkinson shift needed for the
// symqr algorithm. It overwrites the tridiagonal matrix T with Z(T)TZ where
// Z is a product of givens rotations, and returns Z.
Matrix implicitshift(Matrix& T, double PRECISION)
{
  int n = T.nrows(); // It's square
  Matrix Z(n, n);
  double d = (T(n-2, n-2) - T(n-1, n-1))/2.0;
  double u = (d < 0.0 ? -1.0 : 1.0); // Get the sign of d
  u = u*sqrt(d*d + T(n-1, n-2)*T(n-1, n-2));
  u = T(n-1, n-1) -  (T(n-1, n-2)*T(n-1, n-2))/(d + u);
  double x = T(0, 0) - u;
  double z = T(1, 0);
  // Begin main loop
  Vector g; // To store a temporary givens rotation in
  for (int k = 0; k < n-1; k++){
    g = givens(x, z, PRECISION);
    if (k == 0) { // Explicitly get the first givens matrix
      Z = explicitg(g, 0, 1, n);
    } else { // Implicitly keep calculating Z
      Z = givens(Z, g, k, k+1);
    }
    // Calculate TG
    T = givens(T, g, k, k+1);
    // Calculate G(T)T
    T= givens(g, T, k, k+1);
    if (k < n-2){
      x = T(k+1, k);
      z = T(k+2, k);
    }
  }
  return Z;
}


