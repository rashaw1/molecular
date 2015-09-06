/*
 *
 *   PURPOSE: To implement class IntegralEngine, which calculates, stores,
 *            and processes the molecular integrals needed in ab initio
 *            quantum chemistry calculation.
 *
 *   DATE         AUTHOR           CHANGES
 *   ======================================================================
 *   02/09/15     Robert Shaw      Original code.
 *   03/09/15     Robert Shaw      Merged overlap and kinetic integrals.
 *   04/09/15     Robert Shaw      Now supports general contracted.
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
  int M = 0; // No. of spherical basis functions
  for (int i = 0; i < natoms; i++){
    N += m.getAtom(i).getNbfs();
    M += m.getAtom(i).getNSpherical();
  }
  // Cartesian is easy - there are (N^2+N)/2
  // unique 1e integrals and ([(N^2+N)/2]^2 + (N^2+N)/2)/2
  // unique 2e integrals
  int ones = (N*(N+1))/2;
  sizes.resize(4);
  sizes[0] = ones;
  sizes[1] = (ones*(ones+1))/2;
  
  ones = (M*(M+1))/2;
  sizes[2] = ones;
  sizes[3] = (ones*(ones+1))/2;

  formOverlapKinetic();
  formNucAttract();
  std::cout << "\n\nOVERLAP: \n";
  sints.print();
  std::cout << "\n\nKINETIC: \n";
  tints.print();
  std::cout << "\n\n\n";
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

// Sphericalise a shell-pair of 1e- integrals (ints)
// where the shells have angular momentum l1, l2 respectively
// Returns matrix of integrals in canonical order
Matrix IntegralEngine::makeSpherical(int l1, int l2, Matrix& ints) const
{
  // Work out how many m-quantum numbers there are for each l
  int m1n = 2*l1+1;
  int m2n = 2*l2+1;
  
  Matrix retInts(m1n, m2n, 0.0); // The return matrix

  // Calculate all normalisation constants 
  Vector l1norms(m1n); Vector l2norms(m2n);  
  for (int i = -l1; i < l1+1; i++)
    l1norms[i+l1] = sphernorm(l1, i);
  for (int i = -l2; i < l2+1; i++)
    l2norms[i+l2] = sphernorm(l2, i);
  
  // Work out the bounds on summations
  double v1m, v2m, c1, c2;
  int t1max, t2max, v1max, v2max, m1abs, m2abs;

  // Loop over all m1, m2
  for (int m1 = -l1; m1 < l1+1; m1++){
    v1m = (m1 < 0 ? 0.5 : 0.0);
    m1abs = std::abs(m1);
    t1max = std::floor((l1-m1abs)/2.0);
    v1max = std::floor((m1abs/2.0)-v1m);
  
    for (int m2 = -l2; m2 < l2+1; m2++){
      v2m = (m2 < 0 ? 0.5 : 0.0);
      m2abs = std::abs(m2);
      t2max = std::floor((l2-m2abs)/2.0);
      v2max = std::floor((m2abs/2.0)-v2m);
      for (int t1 = 0; t1 < t1max+1; t1++){
	for (int u1 = 0; u1 < t1+1; u1++){
	  for (int v1 = 0; v1 < v1max+1; v1++){
	    // Get clebsch coefficient
	    c1 = clebsch(l1, m1, t1, u1, v1+v1m);
	    int pos1 = ((2*t1+m1abs)*(2*t1+m1abs+1))/2 + 2*t1+m1abs-2*(u1+v1+v1m);
	    for (int t2 = 0; t2 < t2max+1; t2++){
	      for (int u2 = 0; u2 < t2+1; u2++){
		for (int v2 = 0; v2 < v2max+1; v2++){
		  c2 = clebsch(l2, m2, t2, u2, v2+v2m);
		  
		  // Work out position of cartesian integral in ints
		  int pos2 = ((2*t2+m2abs)*(2*t2+m2abs+1))/2 + 2*t2+m2abs-2*(u2+v2+v2m);

		  // Compute the spherical integral
		  retInts(m1+l1, m2+l2) += c1*c2*ints(pos1, pos2);
		}
	      }
	    }
	  }
	}
      }
      // Normalise
      retInts(m1+l1, m2+l2) *= l1norms(m1+l1)*l2norms(m2+l2);
    }
  }
  return retInts;
}
// Same but for 2e- integrals
double IntegralEngine::makeSpherical(int l1, int m1, int l2, int m2,
				     int l3, int m3, int l4, int m4, Matrix& ints) const
{
  double integral = 0.0;
  return integral;
}

// Form the overlap integral matrix sints using the Obara-Saika recurrence relations
// Algorithm:
//  - generate list of basis functions
//  - for each shell  r
//      - for each shell s>=r
//         - loop over prims (u) on in r
//            -loop over prims (v) in s
//                form Si0x,y,z
//                calculate Sij in each cartesian direction
//                Suv = Sij,x*Sij,y*Sij,z
//                form Ti0x,y,z
//                calculate Tij in each direction
//                Tuv = Tij,x*Sij,y*Sij,z + Sij,x*Tij,y*Sij,z + Sij,x*Sij,y*Tij,z
//            end
//         end
//         contract Suv into Smn, and Tuv into Tmn for each pair of basis functions
//         m,n in r, s, respectively.
//      end
//  end 
void IntegralEngine::formOverlapKinetic()
{
  // Get the number of basis functions
  // and the number of shells
  int natoms = molecule.getNAtoms();
  int N = 0; // No. of cartesian cgbfs
  int M = 0; // No. of spherical cgbfs
  int NS = 0;
  for (int i = 0; i < natoms; i++){
    N += molecule.getAtom(i).getNbfs();
    M += molecule.getAtom(i).getNSpherical();
    NS += molecule.getAtom(i).getNshells();
  }
  
  // Resize sints, tints
  sints.assign(M, M, 0.0);
  tints.assign(M, M, 0.0);
  
  // Store the cartesian integrals
  Matrix tempS(N, N); Matrix tempT(N, N);

  // Form a list of basis functions, the atoms they're on,
  // and the shells they belong to 
  Vector atoms(N); Vector bfs(N); Vector shells(N);
  int k = 0;
  Vector temp; 
  for (int i = 0; i < natoms; i++){
    int nbfs = molecule.getAtom(i).getNbfs();
    temp = molecule.getAtom(i).getShells();
    
    for (int j = 0; j < nbfs; j++){
      atoms[k] = i;
      bfs[k] = j;
      int sum = 0; int shell = 0;
      while(sum < j+1){
	sum += temp(shell);
	if(sum < j+1){
	  shell++;
	}
      }
      shells[k] = shell;
      k++;
    }
  }

  int m = 0; // Keep track of basis function count
  int n = 0;

  // Object placeholders
  Vector mcoords; Vector ncoords; Vector mshells; Vector nshells;
  PBF mpbf; PBF npbf; Atom ma; Atom na; 

  // Loop over shells
  for (int r = 0; r < NS; r++){ // Shells on first atom
    // Get the first atom coords and number of prims in this shell
    ma = molecule.getAtom(atoms(m));
    mcoords = ma.getCoords();
    mshells = ma.getShells();
    int mP = ma.getNShellPrims(shells(m));
    
    n = m;
    for (int s = r; s < NS; s++){ // Shells on second atom
      // Get same for second atom
      na = molecule.getAtom(atoms(n));
      ncoords = na.getCoords();
      nshells = na.getShells();
      int nP = na.getNShellPrims(shells(n));

      // Store the primitive integrals
      Matrix overlapPrims(mP, nP);
      Matrix kineticPrims(mP, nP);

      // Loop over primitives
      for (int u = 0; u < mP; u++){
	mpbf = ma.getShellPrim(shells(m), u);

	for (int v = 0; v < nP; v++){
	  npbf = na.getShellPrim(shells(n), v);

	  // Calculate the overlap and kinetic integrals
	  temp = overlapKinetic(mpbf, npbf, mcoords, ncoords);
	  
	  // Store in prim matrices
	  overlapPrims(u, v) = temp(0);
	  kineticPrims(u, v) = temp(1);

	} // End v-loop over prims
      } // End u-loop over prims
      
      // Now we need to contract all the integrals
      // Get shell sizes
      int msize = mshells(shells(m));
      int nsize = nshells(shells(n));

      Vector mplist; Vector nplist;
      Vector mcoeff; Vector ncoeff;
      // Loop

      for (int i = 0; i < msize; i++){
	// Get prim list for this bf, and contraction coeffs
	mplist = ma.getBF(bfs(m+i)).getPrimList();
	mcoeff = ma.getBF(bfs(m+i)).getCoeffs();

	for (int j = i; j < nsize; j++){
	  nplist = na.getBF(bfs(n+j)).getPrimList();
	  ncoeff = na.getBF(bfs(n+j)).getCoeffs();

	  // Form the vector of appropriate prim integrals
	  Vector overInts(mplist.size()*nplist.size());
	  Vector kinInts(mplist.size()*nplist.size());
	  for (int x = 0; x < mplist.size(); x++){
	    for (int y = 0; y < nplist.size(); y++){
	      overInts[x*nplist.size() + y] = overlapPrims(mplist(x), nplist(y));
	      kinInts[x*nplist.size() + y] = kineticPrims(mplist(x), nplist(y));
	    }
	  }
	
	  // Contract cartesian integrals into the temporary matrix
	  // ordered canonically
       	  tempS(m+i, n+j) = makeContracted(mcoeff, ncoeff, overInts);
	  tempT(m+i, n+j) = makeContracted(mcoeff, ncoeff, kinInts);
	  tempS(n+j, m+i) = tempS(m+i, n+j);
	  tempT(n+j, m+i) = tempT(m+i, n+j);
	}
      }
      
      // Increment basis function counts
      n += nsize;
    } // End s-loop over shells

    // Increment basis function counts
    m += mshells(shells(m));
  } // End r-loop over shells
  
}

// Calculate the overlap and kinetic energy integrals between two primitive
// cartesian gaussian basis functions, given the coordinates of their centres
Vector IntegralEngine::overlapKinetic(const PBF& u, const PBF& v, 
				      const Vector& ucoords, const Vector& vcoords) const
{
  Vector rvals(2); // Vector to return answer in

  // Get exponents, norms, and angular momenta
  int ulx = u.getLx(); int uly = u.getLy(); int ulz = u.getLz();
  int vlx = v.getLx(); int vly = v.getLy(); int vlz = v.getLz();
  double unorm = u.getNorm(); double uexp = u.getExponent();
  double vnorm = v.getNorm(); double vexp = v.getExponent();

  // Get the necessary values from getVals
  Vector vals;
  vals = getVals(uexp, vexp, ucoords, vcoords);

  // Store the overlap intermediates for later use
  Matrix Sijx(ulx+1, vlx+1); Matrix Sijy(uly+1, vly+1); Matrix Sijz(ulz+1, vlz+1);

  // Calculate the S00 values in each direction
  double premult = std::sqrt(M_PI/vals(0)); // sqrt(PI/p)
  Sijx(0, 0) = premult*vals(8); // vals(8-10) are the K values
  Sijy(0, 0) = premult*vals(9);
  Sijz(0, 0) = premult*vals(10);

  // Loop to form Si0 in each cartesian direction

  // Use the Obara-Saika recursion formula:
  // S(i+1)j = XPA*Sij + (1/2p)*(i*S(i-1)j + j*Si(j-1))
  // to calculate Si0. 

  // First calculate XPA, YPA, ZPA, 1/2p
  double XPA = vals(2) - ucoords(0);
  double YPA = vals(3) - ucoords(1);
  double ZPA = vals(4) - ucoords(2);
  double one2p = 1.0/(2.0*vals(0)); // 1/2p                                                                                                    

  // Then loop
  double Snext, Scurr, Slast;
  Slast = 0.0;
  Scurr = Sijx(0, 0);
  for (int i = 1; i < ulx+1; i++){
    Snext = XPA*Scurr + one2p*(i-1)*Slast;
    Sijx(i, 0) = Snext;
    Slast = Scurr; Scurr = Snext;
  }
  Slast = 0.0; Scurr = Sijy(0, 0);
  for (int i = 1; i < uly+1; i++){
    Snext = YPA*Scurr + one2p*(i-1)*Slast;
    Sijy(i, 0) = Snext;
    Slast = Scurr; Scurr = Snext;
  }
  Slast = 0.0; Scurr = Sijz(0, 0);
  for (int i = 1; i < ulz+1; i++){
    Snext = ZPA*Scurr + one2p*(i-1)*Slast;
    Sijz(i, 0) = Snext;
    Slast = Scurr; Scurr = Snext;
  }

  // Next we increment the j using the equivalent recursion formula
  // First, calculate XPB, YPB, ZPB
  double XPB = vals(2) - vcoords(0);
  double YPB = vals(3) - vcoords(1);
  double ZPB = vals(4) - vcoords(2);

  // Get the Si1 before looping, if needed
  if(vlx>0){
    for (int k = 0; k < ulx+1; k++){
      int ktemp = (k > 0 ? k-1 : 0); // Avoid out of bounds errors
      Sijx(k, 1) = XPB*Sijx(k,0) + one2p*k*Sijx(ktemp, 0); 
    
      // Then loop
      for (int j = 2; j < vlx+1; j++)
	Sijx(k, j) = XPB*Sijx(k, j-1) + one2p*(k*Sijx(ktemp, j-1) +
					       (j-1)*Sijx(k, j-2));
    }
  }
  // Repeat for y, z
  if(vly>0){
    for (int k = 0; k < uly+1; k++){
      int ktemp= (k > 0 ? k-1 : 0); // Avoid out of bounds errors                                                                                                                                               
      Sijy(k, 1) = YPB*Sijy(k,0) + one2p*k*Sijy(ktemp, 0);

      // Then loop                                                                                                                                                                  
      for (int j = 2; j < vly+1; j++)
        Sijy(k, j) = YPB*Sijy(k, j-1) + one2p*(k*Sijy(ktemp, j-1) +
                                               (j-1)*Sijy(k, j-2));
    }
  }
  if(vlz>0){
    for (int k = 0; k < ulz+1; k++){
      int ktemp= (k > 0 ? k-1 : 0); // Avoid out of bounds errors                                                                                                                                               
      Sijz(k, 1) = ZPB*Sijz(k,0) + one2p*k*Sijz(ktemp, 0);

      // Then loop                                                                                                                                                                             
      for (int j = 2; j < vlz+1; j++)
        Sijz(k, j) = ZPB*Sijz(k, j-1) + one2p*(k*Sijz(ktemp, j-1) +
                                               (j-1)*Sijz(k, j-2));
    }
  }

  // Get final overlap integral
  rvals[0] = unorm*vnorm*Sijx(ulx, vlx)*Sijy(uly, vly)*Sijz(ulz, vlz);
  
  // Now compute the kinetic energy integral
  
  // Start by calculating T00 in each direction
  Matrix Tijx(ulx+1, vlx+1); Matrix Tijy(uly+1, vly+1); Matrix Tijz(ulz+1, vlz+1);
  Tijx(0, 0) = (uexp - 2*uexp*uexp*(XPA*XPA + one2p))*Sijx(0, 0);
  Tijy(0, 0) = (uexp - 2*uexp*uexp*(YPA*YPA + one2p))*Sijy(0, 0);
  Tijz(0, 0) = (uexp - 2*uexp*uexp*(ZPA*ZPA + one2p))*Sijz(0, 0);

  // A couple of repeatedly used multipliers
  double vp = vexp/vals(0); // vexp/p
  double up = uexp/vals(0); // uexp/p

  // Form the Ti0 values
  if (ulx > 0) {
    // Get T10 first
    Tijx(1, 0) = XPA*Tijx(0, 0) + vp*2*uexp*Sijx(1, 0);

    // Loop for rest
    for (int i = 2; i < ulx+1; i++){
      Tijx(i, 0) = XPA*Tijx(i-1, 0) + one2p*(i-1)*Tijx(i-2, 0) + 
	vp*(2*uexp*Sijx(i, 0) - (i-1)*Sijx(i-2, 0));
    }
  }
  // Repeat for y and z components
  if (uly > 0) {
    // Get T10 first
    Tijy(1, 0) = YPA*Tijy(0, 0) + vp*2*uexp*Sijy(1, 0);

    // Loop for rest
    for(int i = 2; i <uly+1; i++){
      Tijy(i, 0) = YPA*Tijy(i-1, 0) + one2p*(i-1)*Tijy(i-2, 0) +
        vp*(2*uexp*Sijy(i, 0) - (i-1)*Sijy(i-2, 0));
    }
  }
  if (ulz > 0) {
    // Get T10 first
    Tijz(1, 0) = ZPA*Tijz(0, 0) + vp*2*uexp*Sijz(1, 0);

    // Loop for rest
    for(int i = 2; i <ulz+1; i++){
      Tijz(i, 0) = ZPA*Tijz(i-1, 0) + one2p*(i-1)*Tijz(i-2, 0) +
        vp*(2*uexp*Sijz(i, 0) - (i-1)*Sijz(i-2, 0));
    }
  }
  
  // Now increment j

  if (vlx > 0){
    for (int k = 0; k < ulx+1; k++){
      int ktemp = (k > 0 ? k-1 : 0);
      Tijx(k, 1) = XPB*Tijx(k,0) + one2p*k*Tijx(ktemp, 0)
	+ up*2*vexp*Sijx(k, 1);
      
      for (int j = 2; j < vlx+1; j++){
	Tijx(k, j) = XPB*Tijx(k, j-1) + one2p*(k*Tijx(ktemp, j-1) + (j-1)*Tijx(k, j-2))
	  + up*(2*vexp*Sijx(k, j) - (j-1)*Sijx(k, j-2));
      }
    }
  }
  // Repeat for y and z
  if (vly > 0){
    for(int k = 0; k <uly+1; k++){
      int ktemp = (k > 0 ? k-1 : 0);
      Tijy(k, 1) = YPB*Tijy(k,0) + one2p*k*Tijy(ktemp, 0)
	+ up*2*vexp*Sijy(k, 1);
      
      for (int j = 2; j< vly+1; j++){
	Tijy(k,j) = YPB*Tijy(k, j-1) +one2p*(k*Tijy(ktemp, j-1) + (j-1)*Tijy(k, j-2))
          + up*(2*vexp*Sijy(k, j) - (j-1)*Sijy(k, j-2));
      } 
    }
  }
  if (vlz > 0){
    for(int k = 0; k <ulz+1; k++){
      int ktemp = (k > 0 ? k-1 : 0);
      Tijz(k, 1) = ZPB*Tijz(k,0) + one2p*k*Tijz(ktemp, 0)
	+ up*2*vexp*Sijz(k, 1);
      
      for (int j = 2; j< vlz+1; j++){
	Tijz(k,j) = ZPB*Tijz(k, j-1) +one2p*(k*Tijz(ktemp, j-1) + (j-1)*Tijz(k, j-2))
          + up*(2*vexp*Sijz(k, j) - (j-1)*Sijz(k, j-2));
      } 
    }
  }

  // Construct the final kinetic energy integral as:
  // Tuv = Tijx*Sijy*Sijz + Sijx*Tijy*Sijz + Sijx*Sijy*Tijz
  rvals[1] = unorm*vnorm*(Tijx(ulx, vlx)*Sijy(uly, vly)*Sijz(ulz, vlz)
    + Sijx(ulx, vlx)*Tijy(uly, vly)*Sijz(ulz, vlz) 
			 + Sijx(ulx, vlx)*Sijy(uly, vly)*Tijz(ulz, vlz));
  
  return rvals;
}

// Form the matrix of nuclear attraction integrals
void IntegralEngine::formNucAttract()
{
}

// Calculate a multipole integral between two bfs a,b 
// about the point c, to a given set of powers in the
// cartesian coordinates of c.
// Uses Obara-Saika recurrence relations. 
double IntegralEngine::multipole(BF& a, BF& b, const Vector& acoords,
				 const Vector& bcoords, const Vector& ccoords,
				 const Vector& powers) const
{
  double integral; // To return answer in
  
  // Get the number of primitives on each, and contraction coefficients
  Vector acoeffs; Vector bcoeffs;
  acoeffs = a.getCoeffs();
  bcoeffs = b.getCoeffs();
  int aN = a.getNPrims();
  int bN = b.getNPrims(); 
  
  // Need to store primitive integrals
  Vector prims(aN*bN);

  // Loop over primitives
  PBF apbf; PBF bpbf;
  for (int u = 0; u < aN; u++){
    apbf = a.getPBF(u);
    for (int v = 0; v < bN; v++){
      bpbf = b.getPBF(v);

      // Calculate the primitive integral
      prims[u*bN+v] = multipole(apbf, bpbf, acoords, bcoords, ccoords, powers);
    }
  }
  
  // Contract
  integral = makeContracted(acoeffs, bcoeffs, prims);
  return integral;
}

// Calculate the above multipole integral between two primitives
double IntegralEngine::multipole(PBF& u,  PBF& v, const Vector& ucoords,
		 const Vector& vcoords, const Vector& ccoords, 
		 const Vector& powers) const
{
  // To be written
  double integral;
  return integral;
}
