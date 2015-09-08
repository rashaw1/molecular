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
 *   06/09/15     Robert Shaw      Nuclear attraction ints for prims.
 *   07/09/15     Robert Shaw      formNucAttract() now works, as does
 *                                 makeSpherical(ints, lnums)
 *   08/09/15     Robert Shaw      Auxiliary 2e- integrals, twoe
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
  std::cout << "OVERLAP: \n"; sints.print(); std::cout << "\n\n";
  formNucAttract();
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

// Sphericalise a matrix of 1e- integrals (ints)
// where the cols have angular momenta lnums.
// Returns matrix of integrals in canonical order
Matrix IntegralEngine::makeSpherical(const Matrix& ints, const Vector& lnums) const
{
  // Calculate the size of matrix needed
  int scount = 0, pcount = 0, dcount = 0, fcount = 0; 
  for (int i = 0; i < lnums.size(); i++){
    switch((int)(lnums(i))){
    case 1: { pcount++; break; } // p-type
    case 2: { dcount++; break; } // d-type
    case 3: { fcount++; break; } // f-type
    default: { scount++; } // Assume s-type
    }
  }

  // Number of spherical basis functions
  int M = scount + pcount + 5*(dcount/6) + 7*(fcount/10); 
  // Number of cartesian basis functions.
  int N = lnums.size();

  // Declare the matrix to return the transformed integrals in
  Matrix retInts;

  // Construct a reduced list of lnums, and
  // corresponding m-quantum numbers
  Vector slnums(M); Vector smnums(M);
  int j=0, k = 0; // Counter for slnums
  while(j < N){
    switch((int)(lnums(j))){
    case 1: { // p-type 
      slnums[k] = 1; smnums[k] = 0;
      slnums[++k] = 1; smnums[k] = -1;
      slnums[++k] = 1; smnums[k++] = 1; 
      j += 3;
      break;
    } 
    case 2: { // d-type
      slnums[k] = 2; smnums[k] = 0;
      slnums[++k] = 2; smnums[k] = -1;
      slnums[++k] = 2; smnums[k] = 1;
      slnums[++k] = 2; smnums[k] = -2;
      slnums[++k] = 2; smnums[k++] = 2;
      j += 6;
      break;
    }
    case 3: { // f-type
      slnums[k] = 3; smnums[k] = 0;
      slnums[++k] = 3; smnums[k] = -1;
      slnums[++k] = 3; smnums[k] = 1;
      slnums[++k] = 3; smnums[k] = -2;
      slnums[++k] = 3; smnums[k] = 2;
      slnums[++k] = 3; smnums[k] = -3;
      slnums[++k] = 3; smnums[k++] = 3;
      j+=10;
      break;
    }
    default: { // assume s-type
      slnums[k] = 0; smnums[k++] = 0;
      j++;
    }
    }
  }

  // Make the transformation matrix
  Matrix trans(M, N, 0.0);
  // Loop through all the slnums, looking up the coefficients
  // as needed. 
  j = 0; // Column and row counters
  for(int i = 0; i < M; i++){
    // Get the coefficients
    formTransMat(trans, i, j, (int)(slnums(i)), (int)(smnums(i)));
    if (slnums(i) - smnums(i) == 0){ // Increment j by a suitable amount
      switch((int)(slnums(i))){
      case 1: { j+=3; break;}
      case 2: { j+=6; break;}
      case 3: { j+=10; break;}
      default: { j+=1; }
      }
    }
  }

  // Now transform the integral matrix
  retInts = trans*ints;
  retInts = retInts*(trans.transpose());
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
//  Transform the integrals to the spherical harmonic basis.
void IntegralEngine::formOverlapKinetic()
{
  // Get the number of basis functions
  // and the number of shells
  int natoms = molecule.getNAtoms();
  int N = 0; // No. of cartesian cgbfs
  int NS = 0; // No. of shells
  for (int i = 0; i < natoms; i++){
    N += molecule.getAtom(i).getNbfs();
    NS += molecule.getAtom(i).getNshells();
  }

  // Resize sints, tints
  sints.assign(N, N, 0.0); tints.assign(N, N, 0.0);
  
  // Form a list of basis functions, the atoms they're on,
  // and the shells they belong to, and the lnums of said bfs 
  Vector atoms(N); Vector bfs(N); Vector shells(N); 
  Vector lnums(N);
  int k = 0;
  Vector temp; Vector temp2;
  for (int i = 0; i < natoms; i++){
    int nbfs = molecule.getAtom(i).getNbfs();
    temp = molecule.getAtom(i).getShells();
    temp2 = molecule.getAtom(i).getLnums();
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
      lnums[k] = temp2(shell);
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

	  // Contract cartesian integrals into the integral matrices
	  // ordered canonically
       	  sints(m+i, n+j) = makeContracted(mcoeff, ncoeff, overInts);
	  tints(m+i, n+j) = makeContracted(mcoeff, ncoeff, kinInts);
	  sints(n+j, m+i) = sints(m+i, n+j);
	  tints(n+j, m+i) = tints(m+i, n+j);
	}
      }
      
      // Increment basis function counts
      n += nsize;
    } // End s-loop over shells

    // Increment basis function counts
    m += mshells(shells(m));
  } // End r-loop over shells
  
  // Transform the matrices to the spherical harmonic basis
  sints = makeSpherical(sints, lnums);
  tints = makeSpherical(tints, lnums);
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
  double integral = 0.0;
  return integral;
}

// Form the matrix of nuclear attraction integrals
// Algorithm:
//    For each pair of atoms m, n:
//       loop over shells, r,  on m
//          loop over shells, s,  on n
//             loop over atomic centres, C
//                loop over primitives u in r
//                   loop over primitivs v in s
//                       calculate the nucAttract between u, v and centre C
//                   end v loop
//                end u loop
//                contract primitive integrals for centre C
//                for each bf pair in (r,s)
//             end C loop
//             sum contributions from each atomic centre for each bf pair in (r,s)
//           end s loop
//        end r loop
//    Transform integrals to spherical harmonic basis.
void IntegralEngine::formNucAttract()
{
    // Get the number of basis functions
  // and the number of shells
  int natoms = molecule.getNAtoms();
  int N = 0; // No. of cartesian cgbfs
  int NS = 0; // No. of shells
  for (int i = 0; i < natoms; i++){
    N += molecule.getAtom(i).getNbfs();
    NS += molecule.getAtom(i).getNshells();
  }

  // Resize naints, and assign all elements to zero
  naints.assign(N, N, 0.0); 
  
  // Form a list of basis functions, the atoms they're on,
  // and the shells they belong to, and get lnums 
  Vector atoms(N); Vector bfs(N); Vector shells(N);
  Vector lnums(N);
  int k = 0;
  Vector temp; Vector temp2;
  for (int i = 0; i < natoms; i++){
    int nbfs = molecule.getAtom(i).getNbfs();
    temp = molecule.getAtom(i).getShells();
    temp2 = molecule.getAtom(i).getLnums();

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
      lnums[k] = temp2(shell);
      k++;
    }
  }

  int m = 0; // Keep track of basis function count
  int n = 0;

  // Object placeholders
  Vector mcoords; Vector ncoords; Vector ccoords;
  Vector mshells; Vector nshells;
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
      Matrix prims(mP, nP);
      
      // Get shell sizes                                                                                                                                                           
      int msize = mshells(shells(m));
      int nsize = nshells(shells(n));
      
      // Loop over atomic centres
      for (int c = 0; c < natoms; c++){
	// Get the coordinates and atomic charge for this centre
	ccoords = molecule.getAtom(c).getCoords();
	int Z = molecule.getAtom(c).getCharge(); 

	// Loop over primitives
	for (int u = 0; u < mP; u++){
	  mpbf = ma.getShellPrim(shells(m), u);
	  
	  for (int v = 0; v < nP; v++){
	    npbf = na.getShellPrim(shells(n), v);
	    
	    // Calculate the nuclear attraction integrals
	    prims(u, v) = nucAttract(mpbf, npbf, mcoords, ncoords, ccoords);
	    	    	  
	  } // End v-loop over prims
	} // End u-loop over prims
      
	// Now we need to contract all the integrals
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
	    Vector ints(mplist.size()*nplist.size());
	    for (int x = 0; x < mplist.size(); x++){
	      for (int y = 0; y < nplist.size(); y++){
		ints[x*nplist.size() + y] = prims(mplist(x), nplist(y));
	      }
	    }
	    
	    // Contract cartesian integrals into the nuclear attraction
	    // matrix, weighting by the atomic charge of centre C
	    naints(m+i, n+j) += -1.0*Z*makeContracted(mcoeff, ncoeff, ints);
	    naints(n+j, m+i) = naints(m+i, n+j);
	  }
	} // End contraction loops
      } // End loop over centres
      // Increment basis function counts
      n += nsize;
    } // End s-loop over shells

    // Increment basis function counts
    m += mshells(shells(m));
  } // End r-loop over shells
  
  // Transform the integrals to the spherical harmonic basis
  naints = makeSpherical(naints, lnums);

}

// Calculate the nuclear attraction integral between two gaussian primitives
// and nucleus C, using the Obara-Saika recurrence relations
// Algorithm:
//    To calculate ijklmn we need to get the aux. integral
//          O(k+l+m+n)_ij0000 
//      followed by 
//          O(m+n)_ijkl00
//      then finally
//          O(0)ijklmn
//   To get each of these, we must start recursively from O(i+j+k+l+m+n)_000000
//      and then increment the first index by vertical recursion,
//      followed by horizontal recursion to increment the second index. This is
//      then repeated at each stage.
double IntegralEngine::nucAttract(const PBF& u, const PBF& v, const Vector& ucoords, 
				  const Vector& vcoords, const Vector& ccoords) const
{
  double integral = 0.0; // To return the answer in

  // Get exponents, norms, and angular momenta
  int ulx = u.getLx(); int uly = u.getLy(); int ulz = u.getLz();
  int vlx = v.getLx(); int vly = v.getLy(); int vlz = v.getLz();
  double unorm = u.getNorm(); double uexp = u.getExponent();
  double vnorm = v.getNorm(); double vexp = v.getExponent();

  // Get the necessary values from getVals
  Vector vals;
  vals = getVals(uexp, vexp, ucoords, vcoords);

  // Determine the maximum N needed for the auxiliary integrals
  // in each cartesian direction
  int Nx = ulx + vlx; int Ny = uly + vly; int Nz = ulz + vlz;
  int N = Nx + Ny + Nz;

  // Calculate/retrieve the needed atomic separations
  double XAB = vals(5); double YAB = vals(6); double ZAB = vals(7);
  double XPA = vals(2) - ucoords(0); double XPC = vals(2) - ccoords(0);
  double YPA = vals(3) - ucoords(1); double YPC = vals(3) - ccoords(1);
  double ZPA = vals(4) - ucoords(2); double ZPC = vals(4) - ccoords(2); 

  // Store the auxiliary integrals
  Matrix aux(N+1, Nx+1, 0.0);
    
  // Calculate the 000000 integral for N, 
  // via the boys function
  double p = vals(0); double pionep = 2.0*M_PI/p; double one2p = 1.0/(2.0*p);
  double K = vals(8)*vals(9)*vals(10); 
  double pRPC2 = p*(XPC*XPC + YPC*YPC + ZPC*ZPC);
  Vector boysval; boysval = boys(pRPC2, N, 0);
  for (int n = 0; n < N+1; n++){ 
    aux(n,0) = pionep*K*boysval(n);
  }

  // Increment the first cartesian index  by the recurrence relation:
  // O(n)_i+1, 0 = XPA*O(n)i,0 - XPC*O(n+1)i,0 + (i/2p)*(O(n)_i-1,0 - O(n+1)_i-1,0)
  // for each n from N-1 to N-Nx
  for (int n = N-1; n > -1 ; n--){
    // Calculate first increment
    aux(n, 1) = XPA*aux(n, 0) - XPC*aux(n+1, 0); 
    int Nmax = (N-n+1 > Nx+1 ? Nx+1 : N-n+1);
    for (int k = 2; k < Nmax; k++){
      aux(n, k) = XPA*aux(n, k-1) - XPC*aux(n+1, k-1) + (k-1)*one2p*(aux(n, k-2) - aux(n+1, k-2));
    }
  }

  // Increment the second index by the horizontal recursion relation if necessary
  for (int n = N-Nx; n > -1; n--){
    for (int j = 1; j < vlx+1; j++){
      for (int i = Nx-1; i > -1; i--){
	aux(n, i) = aux(n, i+1) + XAB*aux(n, i);
      }
    }
  }

  Vector temp; temp = aux.colAsVector(ulx);
  aux.assign(N-Nx+1, Ny+1, 0.0); // Resize, and repeat procedure for next two indices
  aux.setCol(0, temp);

  for (int n = N-Nx-1; n > -1; n--){
    // Calculate first increment
    aux(n, 1) = YPA*aux(n, 0) - YPC*aux(n+1, 0);
    int Nmax = (N-Nx-n+1 > Ny+1 ? Ny+1 : N-Nx-n+1);
    for (int k = 2; k < Nmax; k++){
      aux(n, k) = YPA*aux(n, k-1) - YPC*aux(n+1, k-1) + (k-1)*one2p*(aux(n, k-2) - aux(n+1, k-2));
    }
  }

  // Increment the second index by the horizontal recursion relation if necessary
  for (int n = N-Nx-Ny; n > -1; n--){
    for (int j = 1; j < vly+1; j++){
      for (int i = Ny-1; i > -1; i--){
	aux(n, i) = aux(n, i+1) + YAB*aux(n, i);
      }
    }
  }

  temp = aux.colAsVector(uly);
  aux.assign(Nz+1, Nz+1, 0.0);
  aux.setCol(0, temp);

  for (int n = Nz-1; n > -1; n--){
    // Calculate first increment
    aux(n, 1) = ZPA*aux(n, 0) - ZPC*aux(n+1, 0);
    
    for (int k = 2; k < Nz-n+1; k++){
      aux(n, k) = ZPA*aux(n, k-1) - ZPC*aux(n+1, k-1) + (k-1)*one2p*(aux(n, k-2) - aux(n+1, k-2));
    }
  }

  // Increment the second index by the horizontal recursion relation if necessary
  for (int j = 1; j < vlz+1; j++){
    for (int i = Nz-1; i > -1; i--){
      aux(0, i) = aux(0, i+1) + ZAB*aux(0, i);
    }
  }

  // The final integral is now stored in aux(0, ulz)
  integral = unorm*vnorm*aux(0, ulz);
  return integral;
}

// Calculate the two-electron integrals over
// a set of four primitives u,v,w,x. Forms a matrix of
// [u0|w0] integrals using the Obara-Saika vertical
// and electron transfer recurrence relations. 
// Algorithm:
//   - Form the auxiliary integrals [00|00](m) for
//     m in [0, Lu+Lv+Lw+Lx].
//   - Use the vertical relation to form [u0|00] for
//     u in [0, Lu+Lv+Lw+Lx]
//   - Use electron-transfer relation to form [u0|w0]
//     with u in [Lu, Lu+Lx] and w in [Lw, Lw+Lx]
//   Return these as a matrix ready for contraction, 
//   sphericalisation, and then horizontal recurrence.
Matrix IntegralEngine::twoe(const PBF& u, const PBF& v, const PBF& w,
			    const PBF& x, const Vector& ucoords,
			    const Vector& vcoords, const Vector& wcoords,
			    const Vector& xcoords) const
{
  // Extract all the necessary data from the PBFs
  Vector pvals; Vector qvals;
  pvals = getVals(u.getExponent(), v.getExponent(), ucoords, vcoords);
  qvals = getVals(w.getExponent(), x.getExponent(), wcoords, xcoords);
  
  // Unpack, and calculate distances, exponents, and multipliers.
  double p = pvals(0); double q = wvals(0); double alpha = (p*q)/(p+q);
  double XPA = pvals(2) - ucoords(0); double XPQ = pvals(2) - qvals(2);
  double YPA = pvals(3) - ucoords(1); double YPQ = pvals(3) - qvals(3);
  double ZPA = pvals(4) - ucoords(2); double ZPQ = pvals(4) - qvals(4);
  double XAB = pvals(5); double YAB = pvals(6); double ZAB = pvals(7);
  double XCD = qvals(5); double YCD = qvals(6); double ZCD = qvals(7);
  double K = pvals(8)*pvals(9)*pvals(10)*qvals(8)*qvals(9)*qvals(10);
  double zeromult = 2.0*K*M_PI*M_PI*std::sqrt(M_PI/(p+q))/(p*q);
  double RPQ2 = XPQ*XPQ + YPQ*YPQ + ZPQ*ZPQ;
  double ap = alpha/p; double one2p = 1.0/(2.0*p); double one2q = 1.0/(2.0*q);
  double poq = p/q; 

  // Get the angular momenta
  int Lu = u.getLnum(); int Lv = v.getLnum(); int Lw = w.getLnum();
  int Lx = x.getLnum(); int L = Lu+Lv+Lw+Lx;
 
  int Nx = u.getLx() + v.getLx() + w.getLx() + x.getLx();
  int Ny = u.getLy() + v.getLy() + w.getLy() + x.getLy();
  int Nz = u.getLz() + v.getLz() + w.getLz() + x.getLz();

  // Calculate the O(n)0000;0000;0000 integrals for n in [0, L=Lu+Lv+Lw+Lz]
  // These are given by the formula:
  // [00|00](n) = premult*F_n(alpha*RPQ^2), where F_n is the boys function
  // of order n.
  Vector boysvals(L+1); Matrix aux(L+1, Nx+1, 0.0);

  // First calculate all boys function values
  boysvals = boys(alpha*RPQ2, L, 0);
  // Then convert 
  for (int i = 0; i < L+1; i++)
    aux(i, 0) = premult*boysvals(i);

  // Now we need to use the vertical recurrence relation to increment
  // the first index, one cartesian direction at a time
  for (int n = L-1; n > -1; n--){
    // Calculate first increment
    aux(n, 1) = XPA*aux(n, 0) - ap*XPQ*aux(n+1, 0);
    
    // Now increment as high as needed
    int kmax = (L-n+1 > Nx+1 ? Nx+1 : L-n+1);
    for (int k = 2; k < kmax; k++){
      aux(n, k) = XPA*aux(n, k-1) - ap*XPQ*aux(n+1, k-1) +
	(k-1)*one2p*(aux(n, k-2) - ap*aux(n+1, k-2));
    }
  }

  // for the Nx+1 columns truncated at n = L-Nx, increment
  // the y-direction on u
  
  // Store the results
  Matrix newAux(Nz+1, (Nx+1)*(Ny+1));
  Matrix tempMat(L-Nx+1, Ny+1, 0.0); // For each increment
  for (int m = 0; m < Nx+1; m++){
    // Extract the column
    for (int row = 0; row < L-Nx+1; row++)
      tempMat(row, 0) = aux(row, m);

    // Now increment the y-direction up to Ny
    for (int n = L-Nx-1; n > -1; n--){
      // Calculate the first increment
      tempMat(n, 1) = YPA*tempMat(n, 0) - ap*YPQ*tempMat(n+1, 0);

      int kmax = (L-Nx-n+1 > Ny+1 ? Ny+1 : L-Nx-n+1);
      for (int k = 2; k < Nmax; k++){
	tempMat(n, k) = YPA*tempMat(n, k-1) - ap*YPQ*tempMat(n+1, k-1)+
	  (k-1)*one2p*(tempMat(n, k-2) - ap*tempMat(n+1, k-2));
      }
    }
    
    // Copy the needed columns into newAux
    for (int z = 0; z < Nz+1; z++){
      for (int y = 0; y < Ny+1; y++){
	newAux(z, (m*(Ny+1) + y)) = tempMat(z, y);
      }
    }
  }

  // Next increment in the z-direction for each of these columns of newAux
  aux.assign(1, (Nx+1)*(Ny+1)*(Nz+1)); // To store the results
  tempMat.assign(Nz+1, Nz+1, 0.0); // For each iteration
  for (int m = 0; m < (Nx+1)*(Ny+1); m++){
    // Extract the row from newAux
    for (int row = 0; row < Nz+1; row++)
      tempMat(row, 0) = newAux(row, m);

    // Now increment
    for (int n = Nz-1; n > -1; n--){
      // Calculate first increment
      tempMat(n, 1) = ZPA*tempMat(n, 0) - ap*ZPQ*tempMat(n+1, 0);
      
      for (int k = 2; k < Nz-n+1; k++){
	tempMat(n, k) = ZPA*tempMat(n, k-1) - ap*ZPQ*tempMat(n+1, k-1) +
	  (k-1)*one2p*(tempMat(n, k-2) - ap*tempMat(n+1, k-2));
      }
    }
    
    // And copy results into aux
    for(int z = 0; z < Nz+1; z++)
      aux(0, m*(Nx+1)*(Ny+1) + z) = tempMat(0, z); 
  }

  // We now have all the integrals [u0|00] store in aux
  // for u in [0, L]. They are ordered as follows:
  // [000;0|00], [001;0|00], ..., [00Nz;0|00],
  // [010;0|00], ..., [01Nz; 0|00], ...,
  // [0Ny0; 0| 00], ..., [0NyNz;0|00], ...,
  // [100;0|00], ... and so on.
  
  // Next step is then to use the electron-transfer recurrence relation
  // to transfer cartesian powers from u to w (i.e. from elec 1 to elec 2)
  // We need [u0|w0] for u in [Lu, Lu+Lv] and w in [Lw, Lw+Lx]
  
  double vXxXq = -(v.getExponent()*XAB + x.getExponent()*XCD)/q; 
  int Nyz = (Ny+1)*(Nz+1);

  // Increment all of aux in the x-coordinate up and get for wlx to wlx+xlx
  int wlx = w.getLx(); int xlx = x.getLx();
  nextAux.assign(wlx+xlx+1, (Nx+1)*Nyz, 0.0);
  // Copy in zeroth row from aux, and do the first increment
  nextAux.setRow(0, aux.rowAsVec(0));
  // Set zeroth elements
  for (int col = 0; col < Nyz; col++)
    nextAux(1, col) = vXxXq*nextAux(0, col) - poq*nextAux(0, Nyz+col);
  // Now do the rest of first row
  for (int col1 = 1; col1 < Nx; col1++){
    for (int col2 = 0; col2 < Nyz; col2++){
      nextAux(1, col1*Nyz+col2) = vXxXq*nextAux(0, col1*Nyz+col2) + col1*one2q*nextAux(0, (col1-1)*Nyz+col2) -
	poq*nextAux(0, (col1+1)*Nyz+col2);
    }
  }
  // And the rest of the increments
  for (int m = 2; m < wlx+xlx+1; m++){
    // First portion first
    for (int col = 0; col < Nyz; col++)
      nextAux(m, col) = vXxXq*nextAux(m-1, col) + (m-1)*one2q*nextAux(m-2, col) - poq*nextAux(m-1, Nyz+col);
 
    for (int col1 = 1; col1 < Nx-m+1; col1++){
      for(int col2 = 0; col2 < Nyz; col2++){
	nextAux(m, col1*Nyz+col2) = vXxXq*nextAux(m-1, col1*Nyz+col2) + col1*one2q*nextAux(m-1, (col1-1)*Nyz+col2) +
	  (m-1)*one2q*nextAux(m-2, col1*Nyz+col2) - poq*nextAux(m-1, (col1+1)*Nyz+col2);
      }
    }
  }
 
  // We now have a full set of integrals electron-transferred in the x-coordinate.
  // Do the same for the y-coordinate.  
}
