/*
 *   Implementation of mvector.hpp
 *
 *   DATE               AUTHOR                  CHANGES
 *   ============================================================================
 *   14/08/15           Robert Shaw             Original code.
 *   15/08/15           Robert Shaw             Added error handling.
 *   19/08/15           Robert Shaw             Added p-norm and dot product.
 *   20/08/15           Robert Shaw             Added outer product, angle, sorting.
 *   26/08/15           Robert Shaw             Added cross/triple products.
 */
 
 #include "mvector.hpp"
#include "matrix.hpp"
 #include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
 
// Memory clean up function
void Vector::cleanUp()
{
  v.erase(v.begin(), v.end());
}

// Constructors
Vector::Vector(int length)
{
  n = length; // Set length of vector
  v.resize(length);
}


Vector::Vector(int length, const double& a)
{
  n = length; // Set length of vector
  v.resize(length);
  for (int i = 0; i < length; i++)
    v[i] = a;
}


Vector::Vector(int length, const double* a)
{
  n = length; // Set length
  v.resize(length);
  for (int i = 0; i < length; i++)
    v[i] = a[i];
}

// Copy constructor

Vector::Vector(const Vector& u)
{
  n = u.size(); // Get size
  v = u.v;
}

// Destructor

Vector::~Vector()
{

}

// Shaping functions

void Vector::resize(int length) // Resize with loss of values
{
  n = length; // Reset size
  v.resize(length);
}


void Vector::resizeCopy(int length) 
{ 
  std::vector<double> temp;
  temp = v;
  v.resize(length);
  int max = (n < length ? n : length);
  for (int i = 0; i < max; i++)
    v[i] = temp[i];

  n = length;
}


void Vector::assign(int length, const double& a) // Resizes, setting all vals to a
{
  // Do resizing
  resize(length);
  // Set values to a
  for(int i = 0; i < length; i++)
      v[i] = a; 
}

// Swap elements i and j
void Vector::swap(int i, int j)
{
  double temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}


// Overloaded operators

double& Vector::operator[](int i) // Return v[i]
{
  // No bounds checking
  return v[i];
}


double Vector::operator[](int i) const // Return by value
{
  // No bounds checking
  return v[i];
}


double Vector::operator()(int i) const // Return by value
{
  // No bounds checking
  return v[i];
}


Vector& Vector::operator=(const Vector& u)
{
  int newsize = u.size(); // Get the size
  resize(newsize); // Resize the vector
  // Copy in the values from u
  for (int i = 0; i < newsize; i++){
    v[i] = u.v[i];
  } 
  return *this;
}

// Unary operators

Vector Vector::operator+() const
{
  Vector rvec(n); // Make a vector of size n to return
  for(int i = 0; i < n; i++) { //copy in the values
    rvec[i] = v[i];
  }
  return rvec;
}


Vector Vector::operator-() const
{
  Vector rvec(n); 
  for (int i = 0; i<n; i++) {
    rvec[i] = -v[i];
  }
  return rvec;
}

// Binary operators - will work for vectors of diff. sizes, but throw a warning

Vector Vector::operator+(const Vector& u) const
{
  int size = u.size();
  size = ( n < size ? n : size ); // Set size to the smaller of the vector sizes
  // If different, warn the user, but carry on
  if(size != n) { 
    throw( Error("WARNING", "Vectors are different sizes.") );
  }
  Vector rvec(size); // If one vector is bigger, excess is lost
  for (int i = 0; i < size; i++) { // Do the addition
    rvec[i] = v[i]+u.v[i];
  }
  return rvec;
}


Vector Vector::operator-(const Vector& u) const
{
  int size = u.size();
  size = ( n < size ? n : size ); 
  // If different, warn the user, but carry on                           
  if(size != n) { 
    throw(Error("WARNING", "Vectors are different sizes."));
  } 
  Vector rvec(size); 
  for (int i = 0; i < size; i++) { // Subtract 
    rvec[i] = v[i]-u.v[i]; // This is a left-to-right operator
   }
  return rvec;
}

// Intrinsic functions

// Pretty print
void Vector::print(double PRECISION) const 
{
  double val = 0.0; // Temp printing float, to avoid tiny, tiny numbers
  for (int i = 0; i < n; i++){
    if (fabs(v[i]) > PRECISION){
      val = v[i];
    } else { val = 0.0; }
    std::cout << std::setprecision(12) << std::setw(18) << val;
  }
  std::cout << "\n";
}

// Sort the vector into ascending order, using a quicksort algorithm
void Vector::sort()
{
  if (n > 1) { // Don't bother sorting if only one element (or none!)
  int pivot = rand() % n; // Randomly choose a pivot within length of vector
  swap(n-1, pivot); // Swap the last value in vector with pivot value
  int i, k, p; // Indices for main loop
  i = 0; k = 0; p = n-1;
  // Begin main loop
  while(i < p){
    // Check whether smaller than pivot value
    if (v[i] < v[n-1]) { // Move small values to left
      swap(i, k);
      i++; k++; // Increment i and k
    } else if (v[i] - v[n-1] < 1E-14) { // If they are equal
      swap(i, p-1); // Move equal values to the end 
      p--; 
    } else {
      i++; // Do nothing, but increment counter
    }
  }
  // Now k = the position of the first element in the 'bigger than pivot' list
  // and p = the position of the first element in the 'equal to pivot' list
  // n - p is thus the number of elements equal to the pivot

  // Move pivots to the center
  // m is the size of the smaller of the 'equal to' and 'bigger than' lists
  int m = (p-k < n-p ? p-k : n-p);
  // Copy all the 'equal to' pivots so that first element is at k
  if ( m > 0 ) {
    for(int j = 0; j < m; j++){ 
      swap(j+k, n-m+j);
    }
  }
  // Now k is the position of the first element of the 'equal' to list
  // so k-1 is where the 'smaller than' list ends
  // and the 'bigger than list' starts at k + n - p

  // Sort the subvectors recursively
  // Sort the left list, if it exists
  if (k > 0) {
    // Make temp array
    double* temp = new double[k];
    // Copy values from vector
    for (int a = 0; a < k; a++){
      temp[a] = v[a];
    }
    Vector s1(k, temp);
    delete[] temp;
    s1.sort();
    // Copy values back in
    for (int a = 0; a < k; a++){
      v[a] = s1(a);
    }
  }
  // Sort the right list, if it exists
  if (m > 0) {
    // Make temp array                                            
    double* temp = new double[p-k];
    // Copy values from vector                                               
    for(int a = 0; a <p-k; a++){
      temp[a] =v[n-p+k+a];
    }
    Vector s2(p-k, temp);
    delete[] temp;
    s2.sort();
    // Copy values back in
    for (int a = 0; a < p-k; a++){
      v[n-p+k+a] = s2(a);
    }
  }
  }
}

// Return a sorted array without changing the vector
Vector Vector::sorted() const
{
  Vector u(n); // Return vector
  u = *this;
  u.sort();
  return u;
}
// Friend functions
// Inner (dot) product of two vectors
double inner(const Vector& u, const Vector& w)
{
  double rVal = 0.0; // Return value
  // Get lengths of vectors, check they match
  int usize = u.size();
  int wsize = w.size();
  if(usize == wsize){
    // Calculate inner product
    for (int i = 0; i < usize; i++){
      rVal += u(i)*w(i); // Round brackets return by value, not reference
    }
  } else { // Throw error, and return null vector
    throw( Error("VECDOT", "Vectors different sizes.") );
  }
  return rVal;
}

// Get the outer product of two vectors (a matrix), assuming rhs vector 
// implies its transpose (otherwise it wouldn't work!)
Matrix outer(const Vector& u, const Vector& w)
{
  int usize = u.size();
  int wsize = w.size();
  Matrix rmat(usize, wsize); // Matrix that will be returned
  for (int i = 0; i < usize; i++){
    for (int j = 0; j < wsize; j++){
      // Calculate element ij of rmat
      rmat(i, j) = u(i)*w(j);
    }
  }
  return rmat;
}

// Calculate p-norm of vector u.
// Default to 2-norm, p should be greater than or equal to 0, but no check is given.
// p=0 will give the infinity norm as there isn't an appropriate symbol for infinity
// (and a 0-norm would be pointless).
double pnorm(const Vector& u, int p = 2) 
{
  int usize = u.size();
  double rVal = 0.0; // Initialise return value
  // Check if infinity norm
  if (p == 0){
    // Find the maximum element
    for (int i = 0; i < usize; i++){
      rVal = (fabs(u(i)) > rVal ? fabs(u(i)) : rVal);
    }
  } else {
    for (int i = 0; i < usize; i++) {
      // Calculate (p-norm)^p
      rVal += std::pow(fabs(u(i)), p);   
    }
    rVal = std::pow(rVal, 1.0/(double(p)));
  }
  return rVal;
}

// Get the angle between two vectors
double angle(const Vector& u, const Vector& w)
{
  // Get the magnitudes of the vectors
  double unorm = pnorm(u);
  double wnorm = pnorm(w);
  // Get the dot product
  double dprod = inner(u, w);
  // Use the cosine rule
  // but make sure neither is a zero vector
  double rval = 0.0;
  if(dprod > 1E-12){
    rval = std::acos(dprod/(unorm*wnorm));
  }
  return rval;
}

// !!!ONLY FOR 3D VECTORS!!!

// Cross product
Vector cross(const Vector& u, const Vector& w)
{
  Vector r(3);
  r[0] = u(1)*w(2) - w(1)*u(2);
  r[1] = w(0)*u(2) - w(2)*u(0);
  r[2] = u(0)*w(1) - u(1)*w(0);
  return r;
}

// Triple product
double triple(const Vector& u, const Vector& w, const Vector& z)
{
  Vector r(3);
  r = cross(u, w);
  return inner(r, z);
}
