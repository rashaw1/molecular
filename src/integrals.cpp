/*
 *
 *   PURPOSE: To implement class IntegralEngine, which calculates, stores,
 *            and processes the molecular integrals needed in ab initio
 *            quantum chemistry calculation.
 *
 *   DATE         AUTHOR           CHANGES
 *   ======================================================================
 *   02/09/15     Robert Shaw      Original code.
 *
 */

#include "error.hpp"
#include "integrals.hpp"
#include "mathutil.hpp"
#include "basis.hpp"
#include "logger.hpp"
#include <cmath>
#include <iostream>

// Constructor
IntegralEngine::IntegralEngine(Molecule& m) : molecule(m)
{
  // Calculate sizes
  int natoms = molecule.getNAtoms();
  int N = 0; // No. of cartesian basis functions
  for (int i = 0; i < natoms; i++){
    N += m.getAtom(i).getNbfs();
  }
  // Cartesian is easy - there are (N^2+N)/2
  // unique 1e integrals and ([(N^2+N)/2]^2 + (N^2+N)/2)/2
  // unique 2e integrals
  int ones = (N*(N+1))/2;
  sizes.resize(4);
  sizes[0] = ones;
  sizes[1] = (ones*(ones+1))/2;
  
  // Spherical is much harder - tbc
  sizes[2] = ones;
  sizes[3] = sizes(1);

  formOverlap();
  formKinetic();
}

// Accessors

// Return estimates of the memory that will be needed by the 
// one and two electron integrals. Returns as:
// [1e cart, 2e cart, 1e spher, 2e spher]
Vector IntegralEngine::getEstimates() const
{
  Vector estimates(4);
  // The amount of memory is roughly the number of integrals times the size
  // of a double in memory
  estimates[0] = sizeof(double)*sizes(0);
  estimates[1] = sizeof(double)*sizes(1);
  estimates[2] = sizeof(double)*sizes(2);
  estimates[3] = sizeof(double)*sizes(3);
  return estimates;
}

// Make a matrix of multipole integrals of order pole
Matrix IntegralEngine::makeMultipole(int pole) const
{
  Matrix m;
  return m;
}

// Make a matrix of nuclear attraction integrals
Matrix IntegralEngine::makeNucAttract() const
{
  Matrix na;
  return na;
}

// Make a matrix of electron-electron repulsion integrals
Matrix IntegralEngine::makeERI() const
{
  Matrix eri;
  return eri;
} 

// Utility functions needed to calculate integrals

// Calculates the centre-of-charge coordinates(P), total(p) and reduced(u) exponents, 
// relative coordinates(X, Y, Z), and pre-exponential factors(K) between basis functions
// with exponents a and b, and centres A, B. Vector returned contains:
// (p, u, Px, Py, Pz, X, Y, Z, Kx, Ky, Kz) 
Vector IntegralEngine::getVals(double a, double b, const Vector& A, const Vector& B) const
{
  Vector vals(11); // Return vector
  
  // Calculate p and u
  double p = a+b; 
  double u = (a*b)/(a+b);
  vals[0] = p;
  vals[1] = u;
  
  // Calculate the Ps, XYZ, and Ks
  for (int i = 0; i < 3; i++){
    vals[i+2] = (a*A(i) + b*B(i)) / p; // P
    vals[i+5] = A(i) - B(i); // X, Y, Z
    vals[i+8] = std::exp(-1.0*u*vals(i+5)*vals(i+5)); // K
  }
  
  return vals;
}

// Calculate the spherical normalisation constant for angular and magnetic
// quantum numbers l, m. See Helgaker, Jorgensen, Olsen, Molecular Electronic
// Structure Theory, Chapter 9 pg 338 for formulae.
double IntegralEngine::getN(int l, int m) const
{
  int mabs = std::abs(m);
  double lfact = (double)(fact(l));
  double lplusmfact = (double)(fact(l+mabs));
  double lminmfact = (double)(fact(l-mabs));
  double zerom = (m == 0 ? 2.0 : 1.0);
  
  double N = 1.0/(std::pow(2.0, mabs)*lfact);
  N = N*std::sqrt((2.0*lplusmfact*lminmfact)/zerom);
  return N;
}

// Similarly, get the Clebsch-Gordon coefficient
double IntegralEngine::getC(int l, int m, int t, int u, double v) const
{
  int mabs = std::abs(m);
  double vm = (m < 0 ? 0.5 : 0.0);
  double premult = std::pow(-1.0, t+v-vm) * std::pow(4.0, t);
  double blt = (double)(binom(l, t));
  double blmt = (double)(binom(l-t, mabs+t));
  double btu = (double)(binom(t, u));
  double bm2v = (double)(binom(mabs, 2*v));
  
  return premult*blt*blmt*btu*bm2v;
}

// Contract a set of 1e- integrals
// Assumes that integrals are ordered as: 00, 01, 02, ..., 10, 11, 12, ...,
// and so on, where the first number refers to the index of c1, and the second, 
// that of c2.
double IntegralEngine::makeContracted(Vector& c1, Vector& c2, Vector& ints) const
{
  double integral = 0.0;
  int N1 = c1.size();
  int N2 = c2.size();
  // Loop over contraction coefficients
  for (int i = 0; i < N1; i++){
    for (int j = 0; j < N2; j++){
      integral += c1(i)*c2(j)*ints(i*N2+j);
    }
  }
  return integral;
}

// Do the same but for 2e- integrals
// Assumes integrals are stored as:
// 0000 0001 0002 ... 0010 0011 ....
// 0100 0101 0102 ... 1010 1011 ....
// 0200
// .
// .
// .
// 1000 and so on
double IntegralEngine::makeContracted(Vector& c1, Vector& c2, Vector& c3,
				      Vector& c4, Matrix& ints) const
{
  double integral = 0.0;
  int N1 = c1.size();
  int N2 = c2.size();
  int N3 = c3.size();
  int N4 = c4.size();
  
  // Loop over contraction coefficients
  for (int i = 0; i < N1; i++){
    for (int j = 0; j < N2; j++){
      for (int k = 0; k < N3; k++){
	for (int l = 0; l < N4; l++){
	  integral += c1(i)*c2(j)*c3(k)*c4(l)*ints(i*N2+j, k*N4+l);
	}
      }
    }
  }
  return integral;
}

// Sphericalise a set of 1e- integrals
double IntegralEngine::makeSpherical(int l1, int m1, int l2, int m2, Matrix& ints) const
{
  double integral = 0.0;

  // Get the normalisation 
  double N = getN(l1, m1)*getN(l2, m2);

  // Get the summation limits
  int m1abs = std::abs(m1);
  int m2abs = std::abs(m2);
  double vm1 = (m1 < 0 ? 0.5 : 0);
  double vm2 = (m2 < 0 ? 0.5 : 0);
  double v1lim = std::floor(((double)(m1abs))/2.0 - vm1) + vm1;
  double v2lim = std::floor(((double)(m2abs))/2.0 - vm2) + vm2;
  int t1lim = std::floor((double)(l1-m1abs)/2.0);
  int t2lim = std::floor((double)(l2-m2abs)/2.0);
  
  // Loop over first set of indices
  double C;
  for (int t = 0; t <= t1lim; t++){
    for (int u = 0; u <= t; u++){
      double v = vm1;
      while ( vm1 <= v1lim ){
      }
    }
  }
  return integral;
}

// Same but for 2e- integrals
double IntegralEngine::makeSpherical(int l1, int m1, int l2, int m2,
				     int l3, int m3, int l4, int m4, Matrix& ints) const
{
  double integral = 0.0;
  return integral;
}

// Form the overlap integral matrix sints using the Obara-Saika recurrence relations
// Algorithm: for each pair of atoms A, B
// - shell r on A:
//   - shell s on B:
//     - cgbf a in r:
//       -cgbf b in s:
//         -primitive u in a:
//           -primitive v in b:
//              getVals
//              calculate S00_i = sqrt(PI/p)*Ki in each direction i = x, y, z
//              use S01 and S10, then S11, S21, S12, S22, ... until at correct level
//                with the recurrence relations
//              Suv = product of xyz components
//           end
//         end
//         form vector of primitive integrals, contract
//       end
//     end
//   end
// end
void IntegralEngine::formOverlap()
{
  // Get the number of basis functions
  int natoms = molecule.getNAtoms();
  int N = 0; // No. of cgbfs
  for (int i = 0; i < natoms; i++){
    N += molecule.getAtom(i).getNbfs();
  }
  sints.resize(N, N); // Resize the matrix
  std::cout << "sints is a " << N << " x " << N << " matrix\n";
  
  // Loop over atoms
  Atom ma; Atom na;
  for (int m = 0; m < natoms; m++){
    ma = molecule.getAtom(m);
    // Get the shells and coords for m
    Vector mshells; Vector mcoords;
    mshells = ma.getShells();
    mcoords = ma.getCoords();

    for (int n = m; n < natoms; n++){
      na = molecule.getAtom(n);
      // Get the shells and coords for n
      Vector nshells; Vector ncoords;
      nshells = na.getShells();
      ncoords = na.getCoords();

      // Loop over shells
      int mtotal = 0; // Keep track of how many basis functions have been done
      for (int r = 0; r < mshells.size(); r++){
	BF mbf; BF nbf;
	for (int a = 0; a < mshells(r); a++){

	  mbf = ma.getBF(a+mtotal);
	  // get angular momentum components                                                                                                                                          
	  int alx = mbf.getLx();
	  int aly = mbf.getLy();
	  int alz = mbf.getLz();

	  // Get size of primitive shell                                                                                                                                              
	  int aN = mbf.getNPrims();

	  // Get contraction coefficients                                                                                                                                             
	  Vector acoeffs;
	  acoeffs = mbf.getCoeffs();
	  // Loop over cgbfs
	  int ntotal = 0;
	  for (int s = 0; s < nshells.size(); s++){
	    std::cout << nshells(s) << "\n";
	    for (int b = 0; b < nshells(s); b++){	      

	      nbf = na.getBF(b+ntotal);
	      // Get the components of the angular momenta
	      int blx = nbf.getLx();
	      int bly = nbf.getLy();
	      int blz = nbf.getLz();

	      // Get contraction coefficients
	      Vector bcoeffs; 
	      bcoeffs = nbf.getCoeffs();

	      // Get sizes of primitive shells
     	      int bN = nbf.getNPrims();

	      // Initialise vector to store primitive overlap integrals temporarily
	      Vector prims(aN*bN);

	      // Loop over primitives
	      for (int u = 0; u < aN; u++){
		double uexp = mbf.getPBF(u).getExponent();
		double unorm = mbf.getPBF(u).getNorm();
		
		for (int v = 0; v < bN; v++){
		  double vexp = nbf.getPBF(v).getExponent();
		  double vnorm = nbf.getPBF(v).getNorm();
		  // getVals
		  Vector vals;
		  vals = getVals(uexp, vexp, mcoords, ncoords);
		  
		  // Calculate the S00
		  double premult = std::sqrt(M_PI/vals(0)); // sqrt(PI/p) 
		  Vector Si0x(alx+blx+1);
		  Vector Si0y(aly+bly+1);
		  Vector Si0z(alz+blz+1);
		  Si0x[0] = premult*vals(8); // vals(8-10) are the K values
		  Si0y[0] = premult*vals(9);
		  Si0z[0] = premult*vals(10);
		  
		  // Use the Obara-Saika recursion formula:
		  // S(i+1)j = XPA*Sij + (1/2p)*(i*S(i-1)j + j*Si(j-1))
		  // to calculate Si0. We need to loop up to i = alx+blx, 
		  // aly+bly, alz+blz, to then use horizontal relation
		  // entirely in terms of Si0 members
		  
		  // First calculate XPA, YPA, ZPA, 1/2p
		  double XPA = vals(2) - mcoords(0);
		  double YPA = vals(3) - mcoords(1);
		  double ZPA = vals(4) - mcoords(2);
		  double one2p = 1.0/(2.0*vals(0)); // 1/2p

		  // Then loop
		  double Snext, Scurr, Slast;
		  Slast = 0.0;
		  Scurr = Si0x(0);
		  for (int i = 1; i < alx+2; i++){
		    Snext = XPA*Scurr + one2p*(i-1)*Slast;
		    Si0x[i] = Snext;
		    Slast = Scurr; Scurr = Snext;
		  }
     		  Slast = 0.0; Scurr = Si0y(0);
		  for (int i = 1; i < aly+2; i++){
		    Snext = YPA*Scurr + one2p*(i-1)*Slast;
		    Si0y[i] = Snext;
		    Slast = Scurr; Scurr = Snext;
		  }
		  Slast = 0.0; Scurr = Si0z(0);
		  for (int i = 1; i < alz+2; i++){
		    Snext = ZPA*Scurr + one2p*(i-1)*Slast;
		    Si0z[i] = Snext;
		    Slast = Scurr; Scurr = Snext;
		  }

		  // Next, use the horizontal recurrence relation
		  // Si(j+1) = Si+1,j + XAB*Sij
		  // until blx, bly, blz are reached
		  // This version uses explicitly calculated formulae for Sij in terms
		  // of Si0, and currently only goes up to i = j = 3 (i.e. f-type)
		  // It could easily be extended if there were a need to do so.
		  double Sx, Sy, Sz; // Store the final values of each cartesian component
		  Sx = getS(alx, blx, Si0x, vals(5)); // vals(5-7) are XAB, YAB, ZAB, respectively
		  Sy = getS(aly, bly, Si0y, vals(6));
		  Sz = getS(alz, blz, Si0z, vals(7));
		  
		  // Set the primitive matrix element
		  prims[u*bN+v] = unorm*vnorm*Sx*Sy*Sz;
		} // End primitive loop on b
	      } // End primitive loop on a

	      // Contract
	      sints(mtotal+a, ntotal+b) = makeContracted(acoeffs, bcoeffs, prims);
	      
	    } // End cgbf loop on b

	    ntotal += nshells(s);

	  } // End shell loop on b

	} // End cgbf loop on a

	mtotal += mshells(r);

      } // End shell loop on a      
    } // End atom loop b
  } // End atom loop a
  
  sints.print();

}

// Utility routine to calculate the Sij component from Si0
double IntegralEngine::getS(int l1, int l2, Vector& Si0, double AB) const
{
  double S;
  switch(l2){
  case 1: { // p-type
    S = Si0(l1+l2) + AB*Si0(l1+l2-1);
    break;
  }
  case 2: { // d-type
    S = Si0(l1+l2) + 2*AB*Si0(l1+l2-1) + AB*AB*Si0(l1+l2-2);
    break;
  } 
  case 3: { // f-type
    S = Si0(l1+l2) + 3*AB*Si0(l1+l2-1) + 3*AB*AB*Si0(l1+l2-2) + AB*AB*AB*Si0(l1+l2-3);
    break;
  }
  default: // s-type
    S = Si0(l1);
  }
  return S;
}

// Form the kinetic energy integral matrix tints
void IntegralEngine::formKinetic()
{
}

