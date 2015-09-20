/*
 *
 *   PURPOSE: To implement class IntegralEngine, which calculates, stores,
 *            and processes the molecular integrals needed in ab initio
 *            quantum chemistry calculations.
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
 *   09/09/15     Robert Shaw      Shell 2e- integrals, up to (m0|pq)
 *                                 - need to sphericalise, increment, 
 *                                 then sphericalise again.
 *   10/09/15     Robert Shaw      Finished shell twoe. Got rid of 
 *   							   makeContracted and makeSpherical
 *                                 for 2e ints, as absorbed into twoe.
 */

#include "error.hpp"
#include "integrals.hpp"
#include "mathutil.hpp"
#include "basis.hpp"
#include "logger.hpp"
#include "tensor4.hpp"
#include "tensor6.hpp"
#include "tensor7.hpp"
#include "ten4ten6.hpp"
#include "ten4ten4.hpp"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

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
  int ones = (N*(N+1));
  sizes.resize(4);
  sizes[0] = ones;
  sizes[1] = (ones*(ones+1));
  
  ones = (M*(M+1));
  sizes[2] = ones;
  sizes[3] = (ones*(ones+1));

  molecule.getLog().title("INTEGRAL GENERATION");
  
    molecule.getLog().print("Forming the one electron integrals\n");
  formOverlapKinetic();
  formNucAttract();
    molecule.getLog().print("One electron integrals complete\n");
    molecule.getLog().localTime();
    
    Vector ests = getEstimates();
          
    if (molecule.getLog().direct()){
        molecule.getLog().print("Two electron integrals to be calculated on the fly.\n");
    } else if(molecule.getLog().getMemory() > ests(3)){ // Check memory requirements
        formERI(false); // Don't write to file
        if (molecule.getLog().twoprint()) {
            printERI(molecule.getLog().getIntFile(), M);
        }
    } else {
        molecule.getLog().print("Writing the two electron integrals to file\n");
        formERI(true); // Write to file without forming twoints
    }
    //std::cout << "\n\n";
    //sints.print();
    //std::cout << "\n\n";
	//tints.print(); std::cout << "\n";
    //naints.print(); std::cout << "\n";
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
  double TOMB = 1.0/(1024.0 * 1024.0);
  estimates[0] = TOMB*sizeof(double)*sizes(0);
  estimates[1] = TOMB*sizeof(double)*sizes(1);
  estimates[2] = TOMB*sizeof(double)*sizes(2);
  estimates[3] = TOMB*sizeof(double)*sizes(3);
  return estimates;
}

// Return a particular integral from twoints (taking into account symmetries)
double IntegralEngine::getERI(int i, int j, int k, int l) const
{
  double rval = 0.0;
  if (twoints(i, j, k, l) != 0.0 ){ rval = twoints(i, j, k, l); }
  else if (twoints(i, j, l, k) != 0.0){ rval = twoints(i, j, l, k); }
  else if (twoints(j, i, k, l) != 0.0){ rval = twoints(j, i, k, l); }
  else if (twoints(j, i, l, k) != 0.0){ rval = twoints(j, i, l, k); }
  else if (twoints(k, l, i, j) != 0.0){ rval = twoints(k, l, i, j); }
  else if (twoints(k, l, j, i) != 0.0){ rval = twoints(k, l, j, i); }
  else if (twoints(l, k, i, j) != 0.0){ rval = twoints(l, k, i, j); }
  else { rval = twoints(l, k, j, i); }
  return rval;
} 

// Form a tensor of the two-electron integrals (only call if there is
// definitely enough memory!)
void IntegralEngine::formERI(bool tofile)
{
    // Get the number of basis functions
    // and the number of shells
    int natoms = molecule.getNAtoms();
    int N = 0; // No. of cartesian cgbfs
    int NS = 0; // No. of shells
    int NSpher = 0; // No. of spherical cgbfs
    for (int i = 0; i < natoms; i++){
        N += molecule.getAtom(i).getNbfs();
        NS += molecule.getAtom(i).getNshells();
        NSpher += molecule.getAtom(i).getNSpherical();
    }
    
    // Form a list of basis functions, the atoms they're on,
    // and the shells they are in.
    Vector atoms(N); Vector bfs(N); Vector shells(N);
    int k = 0;
    Vector temp; Vector temp2;
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
    
    // Counters to keep track of cart. bfs
    int m = 0; int n = 0; int p = 0; int q = 0;
    // Counters to keep track of spher. bfs
    int a = 0; int b = 0; int c = 0; int d = 0;
    
    Atom ma; Atom na; Atom pa; Atom qa;
    // Loop over shell quartets
    Tensor4 tempInts;
    
    if (!tofile){
        twoints.assign(NSpher, NSpher, NSpher, NSpher, 0.0);
    }
    
    std::ofstream& output = molecule.getLog().getIntFile();
    // Do the "diagonal" elements for pre-screening purposes
    Matrix prescreen(NS, NS, 0.0);
    for (int r = 0; r < NS; r++) {
        ma = molecule.getAtom(atoms(m));
        Vector mshells = ma.getShells();
        int spherR = ma.getNSpherShellBF(shells(m));
        
        n = m;
        b = a;
        for (int s = r; s < NS; s++) {
            na = molecule.getAtom(atoms(n));
            Vector nshells = na.getShells();
            int spherS = na.getNSpherShellBF(shells(n));
            
            // Get the integrals
            tempInts = twoe(ma, na, ma, na, shells(m), shells(n), shells(m), shells(n));
            
            // Copy into two ints, and find maximum element
            double maxval = 0.0;
            double tempval;
            for (int w = 0; w < spherR; w++){
                for (int x = 0; x < spherS; x++){
                    for (int y = 0; y < spherR; y++){
                        for (int z = 0; z < spherS; z++){
                            if (w==y && x==z){
                                tempval = fabs(tempInts(w,x,y,z));
                                maxval = (tempval > maxval ? tempval : maxval);
                            }
                            if (!tofile) {
                                twoints(a+w, b+x, a+y, b+z) = tempInts(w, x, y, z);
                            } else {
                                output << std::setw(6) << a+w+1;
                                output << std::setw(6) << b+x+1;
                                output << std::setw(6) << a+y+1;
                                output << std::setw(6) << b+z+1;
                                output << std::setw(20) << tempInts(w, x, y, z);
                                output << "\n";
                            }
                        }
                    }
                }
            }
            
            // Put maxval in the prescreening matrx
            prescreen(r, s) = std::sqrt(maxval);
            prescreen(s, r) = prescreen(r, s);
            
            // Increment no. of bfs
            n += nshells(shells(n));
            b += spherS;
        }
        
        m += mshells(shells(m));
        a += spherR;
    }
    
    // Update the log file
    molecule.getLog().print("Forming the two electron repulsion integrals.\n");
    molecule.getLog().print("PRESCREENING MATRIX:\n");
    molecule.getLog().print(prescreen);
    molecule.getLog().print("\n\n");
    
    m = n = a = b = 0;
    for (int r = 0; r < NS; r++){
        ma = molecule.getAtom(atoms(m));
        Vector mshells = ma.getShells();
        int spherR = ma.getNSpherShellBF(shells(m));
        
        n = m;
        b = a;
        for (int s = r; s < NS; s++){
            na = molecule.getAtom(atoms(n));
            Vector nshells = na.getShells();
            int spherS = na.getNSpherShellBF(shells(n));
            
            p = m;
            c = a;
            for (int t = r; t < NS; t++){
                pa = molecule.getAtom(atoms(p));
                Vector pshells = pa.getShells();
                int spherT = pa.getNSpherShellBF(shells(p));
                
                q = p;
                d = c;
                for (int u = t; u < NS; u++){
                    qa = molecule.getAtom(atoms(q));
                    Vector qshells = qa.getShells();
                    int spherU = qa.getNSpherShellBF(shells(q));
                    
                    // Prescreen before doing the integrals - Cauchy-Schwarz
                    if ( ((r!=t) || (s!=u)) && (prescreen(r, s)*prescreen(t, u) > molecule.getLog().thrint()) ) {
                        // Get the integrals
                        tempInts = twoe(ma, na, pa, qa, shells(m), shells(n),
                                        shells(p), shells(q));
                        
                        // Copy them into two ints
                        for (int w = 0; w < spherR; w++){
                            for (int x = 0; x < spherS; x++){
                                for (int y = 0; y < spherT; y++){
                                    for (int z = 0; z < spherU; z++){
                                        if (!tofile){
                                            twoints(a+w, b+x, c+y, d+z) = tempInts(w, x, y, z);
                                        } else {
                                            output << std::setw(6) << a+w+1;
                                            output << std::setw(6) << b+x+1;
                                            output << std::setw(6) << c+y+1;
                                            output << std::setw(6) << d+z+1;
                                            output << std::setw(20) << tempInts(w, x, y, z);
                                            output << "\n";
                                        }
                                    } // end z-loop
                                } // end y-loop
                            } // end x-loop
                        } // end w-loop
                    }
                    
                    q += qshells(shells(q));
                    d += spherU;
                } // end u-loop
                
                p += pshells(shells(p));
                c += spherT;
            } // end t-loop
            
            n += nshells(shells(n));
            b += spherS;
        } // end s-loop
        
        m += mshells(shells(m));
        a += spherR;
    } // end r-loop
    
    molecule.getLog().print("Two electron integrals completed.\n");
    if(!tofile){
        std::string mem = "Approximate memory usage = ";
        mem += std::to_string(NSpher*NSpher*NSpher*NSpher*sizeof(double)/(1024.0*1024.0));
        molecule.getLog().print(mem);
    }
    molecule.getLog().localTime();
}

// Print a sorted list of ERIs to ostream output
void IntegralEngine::printERI(std::ostream& output, int NSpher) const
{
  output << " TWO ELECTRON INTEGRALS: \n";
  // print them out
  int icount = 0;
  int scount = 0;
  for (int c1 = 0; c1 < NSpher; c1++){
    for (int c2 = 0; c2 < c1+1; c2++){
      for (int c3 = 0; c3 < c1+1; c3++){
        for(int c4 = 0; c4 < c3+1; c4++){
            icount++;
            double multiplier = 0.125;
            if (c1!=c2) { multiplier *= 2.0; }
            if (c3!=c4) { multiplier *= 2.0; }
            if ((c1+c2) != (c3+c4)) { multiplier*=2.0; }
	    else if ((c1!=c3) && (c1 != c4)) { multiplier *=2.0; }
            else if ((c2!=c3) && (c2 != c4)) { multiplier *=2.0; }
            if (fabs(getERI(c4, c3, c2, c1)) < molecule.getLog().thrint()) { scount++; multiplier = 0; }
            output << std::setw(6) << c1+1;
            output << std::setw(6) << c2+1;
            output << std::setw(6) << c3+1;
            output << std::setw(6) << c4+1;
            output << std::setw(20) << multiplier*getERI(c4, c3, c2, c1);
            output << "\n";
        }
      }
    }
  }			
  output << "N 2e Ints: " << icount << "\n";
  output << "N insig. Ints: " << scount << "\n";
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
      slnums[k] = 1; smnums[k] = 1;
      slnums[++k] = 1; smnums[k] = -1;
      slnums[++k] = 1; smnums[k++] = 0; 
      j += 3;
      break;
    } 
    case 2: { // d-type
      slnums[k] = 2; smnums[k] = 0;
      slnums[++k] = 2; smnums[k] = -2;
      slnums[++k] = 2; smnums[k] = 1;
      slnums[++k] = 2; smnums[k] = 2;
      slnums[++k] = 2; smnums[k++] = -1;
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
  j = 0; // Column counter
  int m = 0; // m-counter
  for(int i = 0; i < M; i++){
    // Get the coefficients
    formTransMat(trans, i, j, (int)(slnums(i)), (int)(smnums(i)));
    if (m == 2*slnums(i)){ // Increment j by a suitable amount
      switch((int)(slnums(i))){
      case 1: { j+=3; break;}
      case 2: { j+=6; break;}
      case 3: { j+=10; break;}
      default: { j+=1; }
      }
      m = 0;
    } else { m++; }
  }

  // Now transform the integral matrix
  retInts = trans*ints;
  retInts = retInts*(trans.transpose());
  return retInts;
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

	for (int j = 0; j < nsize; j++){
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
      Matrix prims(mP, nP, 0.0);
      
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
	  
	  for (int j = 0; j < nsize; j++){
	    nplist = na.getBF(bfs(n+j)).getPrimList();
	    ncoeff = na.getBF(bfs(n+j)).getCoeffs();
	    
	    // Form the vector of appropriate prim integrals
	    Vector ints(mplist.size()*nplist.size(), 0.0);
	    for (int x = 0; x < mplist.size(); x++){
	      for (int y = 0; y < nplist.size(); y++){
		ints[x*nplist.size() + y] = prims(mplist(x), nplist(y));
	      }
	    }
	    // Contract cartesian integrals into the nuclear attraction
	    // matrix, weighting by the atomic charge of centre Cf 	
	    naints(m+i, n+j) += -1.0*Z*makeContracted(mcoeff, ncoeff, ints);
	  }
	} // End contraction loops
      } // End loop over centres
      // Increment basis function counts
      n += nsize;
    } // End s-loop over shells

    // Increment basis function counts
    m += mshells(shells(m));
  } // End r-loop over shells
  
  // Symmetrise
  for (int i = 0; i < N; i++){
    for (int j = i; j < N; j++){
      naints(j, i) = naints(i, j);
    }
  }
  
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
  Tensor7 Aux(N+1, Nx+1, Nx+1, Ny+1, Ny+1, Nz+1, Nz+1, 0.0);
  
  // Calculate the 000000 integral for N, 
  // via the boys function
  double p = vals(0); double pionep = 2.0*M_PI/p; double one2p = 1.0/(2.0*p);
  double K = vals(8)*vals(9)*vals(10); 
  double pRPC2 = p*(XPC*XPC + YPC*YPC + ZPC*ZPC);
  Vector boysval; boysval = boys(pRPC2, N, 0);
  for (int n = 0; n < N+1; n++){ 
  	Aux(n, 0, 0, 0, 0, 0, 0) = pionep*K*boysval(n);
  }

  // Increment the first cartesian index  by the recurrence relation:
  // O(n)_i+1, 0 = XPA*O(n)i,0 - XPC*O(n+1)i,0 + (i/2p)*(O(n)_i-1,0 - O(n+1)_i-1,0)
  // for each n from N-1 to N-Nx
  if ( Nx > 0) {
  for (int n = N-1; n > -1 ; n--){
    // Calculate first increment
    Aux(n, 1, 0, 0, 0, 0, 0) = XPA*Aux(n, 0, 0, 0, 0, 0, 0) - XPC*Aux(n+1, 0, 0, 0, 0, 0, 0);
    int Nmax = (N-n+1 > Nx+1 ? Nx+1 : N-n+1);
    for (int k = 2; k < Nmax; k++){
      Aux(n, k, 0, 0, 0, 0, 0) = XPA*Aux(n, k-1, 0, 0, 0, 0, 0) - XPC*Aux(n+1, k-1, 0, 0, 0, 0, 0)
      	+(k-1)*one2p*(Aux(n, k-2, 0, 0, 0, 0, 0) - Aux(n+1, k-2, 0, 0, 0, 0, 0));
    }
  }
  }
  
  // Increment the second index by the horizontal recursion relation if necessary
  for (int n = N-Nx; n > -1; n--){
    for (int j = 1; j < vlx+1; j++){
      for (int i = Nx-j; i > -1; i--){
		Aux(n, i, j, 0, 0, 0, 0) = Aux(n, i+1, j-1, 0, 0, 0, 0) + XAB*Aux(n, i, j-1, 0, 0, 0, 0);
      }
    }
  }
  
  if( Ny > 0 ) {
  for (int n = N-Nx-1; n > -1; n--){
    // Calculate first increment
    Aux(n, ulx, vlx, 1, 0, 0, 0) = YPA*Aux(n, ulx, vlx, 0, 0, 0, 0) - YPC*Aux(n+1, ulx, vlx, 0, 0, 0, 0);
    int Nmax = (N-Nx-n+1 > Ny+1 ? Ny+1 : N-Nx-n+1);
    for (int k = 2; k < Nmax; k++){
      Aux(n, ulx, vlx, k, 0, 0, 0) = YPA*Aux(n, ulx, vlx, k-1, 0, 0, 0) - YPC*Aux(n+1, ulx, vlx, k-1, 0, 0, 0)
      	+(k-1)*one2p*(Aux(n, ulx, vlx, k-2, 0, 0, 0) - Aux(n+1, ulx, vlx, k-2, 0, 0, 0));
    }
  }
  }
  
  // Increment the second index by the horizontal recursion relation if necessary
  for (int n = N-Nx-Ny; n > -1; n--){
    for (int j = 1; j < vly+1; j++){
      for (int i = Ny-j; i > -1; i--){
      	Aux(n, ulx, vlx, i, j, 0, 0) = Aux(n, ulx, vlx, i+1, j-1, 0, 0) + YAB*Aux(n, ulx, vlx, i, j-1, 0, 0);
      }
    }
  }

  for (int n = Nz-1; n > -1; n--){
    // Calculate first increment
    Aux(n, ulx, vlx, uly, vly, 1, 0) = ZPA*Aux(n, ulx, vlx, uly, vly, 0, 0) - ZPC*Aux(n+1, ulx, vlx, uly, vly, 0, 0);
    
    for (int k = 2; k < Nz-n+1; k++){
      Aux(n, ulx, vlx, uly, vly, k, 0) = ZPA*Aux(n, ulx, vlx, uly, vly, k-1, 0)
      	-ZPC*Aux(n+1, ulx, vlx, uly, vly, k-1, 0) + (k-1)*one2p*(Aux(n, ulx, vlx, uly, vly, k-2, 0) - Aux(n+1, ulx, vlx, uly, vly, k-2, 0));
    }
  }

  // Increment the second index by the horizontal recursion relation if necessary
  for (int j = 1; j < vlz+1; j++){
    for (int i = Nz-j; i > -1; i--){
	  Aux(0, ulx, vlx, uly, vly, i, j) = Aux(0, ulx, vlx, uly, vly, i+1, j-1) + ZAB*Aux(0, ulx, vlx, uly, vly, i, j-1);
    }
  }

  // The final integral is now stored in aux(0, ulz)
  integral = unorm*vnorm*Aux(0, ulx, vlx, uly, vly, ulz, vlz);
  return integral;
}

Tensor4 IntegralEngine::makeE(int u, int v, double K, double p, double PA, double PB) const
{
	int N = u+v;
	Tensor4 E(N+1, N+1, N+1, 1, 0.0);
	
	// Set initial values
	E(0, 0, 0, 0) = K;
	
	double one2p = 1.0/(2.0*p);
	
	if (N > 0) {
		E(0, 1, 0, 0) = PA*E(0, 0, 0, 0);
		E(0, 0, 1, 0) = PB*E(0, 0, 0, 0);
		E(1, 1, 0, 0) = one2p*E(0, 0, 0, 0);
		E(1, 0, 1, 0) = one2p*E(0, 0, 0, 0);
		E(1, 1, 1, 0) = one2p*(E(0, 0, 1, 0) + E(0, 1, 0, 0));
		E(0, 1, 1, 0) = PB*E(0, 1, 0, 0) + E(1, 1, 0, 0);
	}
	
	if (N > 1) {
		E(0, 2, 0, 0) = PA*E(0, 1, 0, 0) + E(1, 1, 0, 0);
		E(0, 0, 2, 0) = PB*E(0, 0, 1, 0) + E(1, 0, 1, 0);
		E(1, 2, 0, 0) = one2p*2*E(0, 1, 0, 0);
		E(1, 0, 2, 0) = one2p*2*E(0, 0, 1, 0);
		E(1, 2, 1, 0) = one2p*(2*E(0, 1, 1, 0) + E(0, 2, 0, 0));
		E(1, 1, 2, 0) = one2p*(E(0, 0, 2, 0) + 2*E(0, 1, 1, 0));
		E(0, 2, 1, 0) = PA*E(0, 1, 1, 0) + E(1, 1, 1, 0);
		E(0, 1, 2, 0) = PB*E(0, 1, 1, 0) + E(1, 1, 1, 0);
		E(0, 2, 2, 0) = PA*E(0, 1, 2, 0) + E(1, 1, 2, 0);
		E(1, 2, 2, 0) = one2p*2*(E(0, 1, 2, 0) + E(0, 2, 1, 0));
		E(2, 2, 0, 0) = one2p*E(1, 1, 0, 0);
		E(2, 0, 2, 0) = one2p*E(1, 0, 1, 0);
		E(2, 1, 1, 0) = (one2p/2.0)*(E(1, 0, 1, 0) + E(1, 1, 0, 0));
		E(2, 2, 1, 0) = (one2p/2.0)*(2*E(1, 1, 1, 0) + E(1, 2, 0, 0));
		E(2, 1, 2, 0) = (one2p/2.0)*(E(1, 0, 2, 0) + 2*E(1, 1, 1, 0));
		E(2, 2, 2, 0) = one2p*(E(1, 1, 2, 0) + E(1, 2, 1, 0));
	}	
	
	if (N > 2) {
		E(0, 3, 0, 0) = PA*E(0, 2, 0, 0) + E(1, 2, 0, 0);
		E(0, 0, 3, 0) = PB*E(0, 0, 2, 0) + E(1, 0, 2, 0);
		E(1, 3, 0, 0) = one2p*2*E(0, 2, 0, 0);
		E(1, 0, 3, 0) = one2p*2*E(0, 0, 2, 0);
		E(1, 3, 1, 0) = one2p*(3*E(0, 2, 1, 0) + E(0, 3, 0, 0));
		E(1, 1, 3, 0) = one2p*(E(0, 0, 3, 0) + E(0, 1, 2, 0));
		E(0, 3, 1, 0) = PA*E(0, 2, 1, 0) + E(1, 2, 1, 0);
		E(0, 1, 3, 0) = PB*E(0, 1, 2, 0) + E(1, 1, 2, 0);
		E(1, 3, 2, 0) = one2p*(3*E(0, 2, 2, 0) + 2*E(0, 3, 1, 0));
		E(1, 2, 3, 0) = one2p*(2*E(0, 1, 3, 0) + 3*E(0, 2, 2, 0));
		E(0, 3, 2, 0) = PA*E(0, 2, 2, 0) + E(1, 2, 2, 0);
		E(0, 2, 3, 0) = PB*E(0, 2, 2, 0) + E(1, 2, 2, 0);
		E(1, 3, 3, 0) = one2p*3*(E(0, 2, 3, 0)  + E(0, 3, 2, 0));
		E(0, 3, 3, 0) = PA*E(0, 2, 3, 0) + E(1, 2, 3, 0);
		E(2, 3, 0, 0) = (one2p/2.0)*3*E(1, 2, 0, 0);
		E(2, 0, 3, 0) = (one2p/2.0)*3*E(1, 0, 2, 0);
		E(2, 3, 1, 0) = (one2p/2.0)*(3*E(1, 2, 1, 0) + E(1, 3, 0, 0));
		E(2, 1, 3, 0) = (one2p/2.0)*(E(1, 0, 3, 0) + 3*E(1, 1, 2, 0));
		E(2, 2, 3, 0) = (one2p/2.0)*(2*E(1, 1, 3, 0) + 3*E(1, 2, 2, 0));
		E(2, 3, 2, 0) = (one2p/2.0)*(3*E(1, 2, 2, 0) + 2*E(1, 3, 1, 0));
		E(2, 3, 3, 0) = (one2p/2.0)*3*(E(1, 2, 3, 0) + E(1, 3, 2, 0));
		E(3, 3, 0, 0) = one2p*E(2, 2, 0, 0);
		E(3, 0, 3, 0) = one2p*E(2, 0, 2, 0);
		E(3, 3, 1, 0) = (one2p/3.0)*(3*E(2, 2, 1, 0) + E(2, 3, 0, 0));
		E(3, 1, 3, 0) = (one2p/3.0)*(E(2, 0, 3, 0) + 3*E(2, 1, 2, 0));
		E(3, 3, 2, 0) = (one2p/3.0)*(3*E(2, 2, 2, 0) + 2*E(2, 3, 1, 0));
		E(3, 2, 3, 0) = (one2p/3.0)*(2*E(2, 1, 3, 0) + 3*E(2, 2, 2, 0));
		E(3, 3, 3, 0) = one2p*(E(2, 2, 3, 0) + E(2, 3, 2, 0));
	}
	
	if (N>3) {
		E(0, 4, 0, 0) = PA*E(0, 3, 0, 0) + E(1, 3, 0, 0);
		E(0, 0, 4, 0) = PB*E(0, 0, 3, 0) + E(1, 0, 3, 0);
		E(1, 4, 0, 0) = one2p*4*E(0, 3, 0, 0);
		E(1, 0, 4, 0) = one2p*4*E(0, 0, 3, 0);
		E(1, 4, 1, 0) = one2p*(4*E(0, 3, 1, 0) + E(0, 4, 0, 0));
		E(1, 1, 4, 0) = one2p*(E(0, 0, 4, 0) + 4*E(0, 1, 3, 0));
		E(0, 4, 1, 0) = PA*E(0, 3, 1, 0) + E(1, 3, 1, 0);
		E(0, 1, 4, 0) = PB*E(0, 1, 3, 0) + E(1, 1, 3, 0);
		E(0, 4, 2, 0) = PA*E(0, 3, 2, 0) + E(1, 3, 2, 0);
		E(0, 2, 4, 0) = PB*E(0, 2, 3, 0) + E(1, 2, 3, 0);
		E(1, 4, 2, 0) = one2p*2*(2*E(0, 3, 2, 0) + E(0, 4, 1, 0));
		E(1, 2, 4, 0) = one2p*2*(E(0, 1, 4, 0) + 2*E(0, 2, 3, 0));
		E(1, 4, 3, 0) = one2p*(4*E(0, 3, 3, 0) + 3*E(0, 4, 2, 0));
		E(1, 3, 4, 0) = one2p*(3*E(0, 2, 4, 0) + 4*E(0, 3, 3, 0));
		E(0, 4, 3, 0) = PA*E(0, 3, 3, 0) + E(1, 3, 3, 0);
		E(0, 3, 4, 0) = PB*E(0, 3, 3, 0) + E(1, 3, 3, 0);
		E(1, 4, 4, 0) = one2p*4*(E(0, 3, 4, 0) + E(0, 4, 3, 0));
		E(0, 4, 4, 0) = PA*E(0, 3, 4, 0) + E(1, 3, 4, 0);
		E(2, 4, 0, 0) = one2p*2*E(1, 3, 0, 0);
		E(2, 0, 4, 0) = one2p*2*E(1, 0, 3, 0);
		E(2, 4, 1, 0) = (one2p/2.0)*(4*E(1, 3, 1, 0) + E(1, 4, 0, 0));
		E(2, 1, 4, 0) = (one2p/2.0)*(E(1, 0, 4, 0) + 4*E(1, 1, 3, 0));
		E(2, 4, 2, 0) = one2p*(2*E(1, 3, 2, 0) + E(1, 4, 1, 0));
		E(2, 2, 4, 0) = one2p*(E(1, 1, 4, 0) + 2*E(1, 2, 3, 0));
		E(2, 4, 3, 0) = (one2p/2.0)*(4*E(1, 3, 3, 0) + 3*E(1, 4, 2, 0));
		E(2, 3, 4, 0) = (one2p/2.0)*(3*E(1, 2, 4, 0) + 4*E(1, 3, 3, 0));
		E(2, 4, 4, 0) = one2p*2*(E(1, 3, 4, 0) + E(1, 4, 3, 0));
		E(3, 4, 0, 0) = (one2p/3.0)*4*E(2, 3, 0, 0);
		E(3, 0, 4, 0) = (one2p/3.0)*4*E(2, 0, 3, 0);
		E(3, 4, 1, 0) = (one2p/3.0)*(4*E(2, 3, 1, 0) + E(2, 4, 0, 0));
		E(3, 1, 4, 0) = (one2p/3.0)*(E(2, 0, 4, 0) + E(2, 1, 3, 0));
		E(3, 4, 2, 0) = (one2p/3.0)*2*(2*E(2, 3, 2, 0) + E(2, 4, 1, 0));
		E(3, 2, 4, 0) = (one2p/3.0)*2*(E(2, 1, 4, 0) + 2*E(2, 2, 3, 0));
		E(3, 4, 3, 0) = (one2p/3.0)*(4*E(2, 3, 3, 0) + 3*E(2, 4, 2, 0));
		E(3, 3, 4, 0) = (one2p/3.0)*(3*E(2, 2, 4, 0) + 4*E(2, 3, 3, 0));
		E(3, 4, 4, 0) = (one2p/3.0)*4*(E(2, 3, 4, 0) + E(2, 4, 3, 0));
		E(4, 4, 0, 0) = one2p*E(3, 3, 0, 0);
		E(4, 0, 4, 0) = one2p*E(3, 0, 3, 0);
		E(4, 4, 1, 0) = (one2p/4.0)*(4*E(3, 3, 1, 0) + E(3, 4, 0, 0));
		E(4, 1, 4, 0) = (one2p/4.0)*(E(3, 0, 4, 0) + 4*E(3, 1, 3, 0));
		E(4, 4, 2, 0) = (one2p/2.0)*(2*E(3, 3, 2, 0) + E(3, 4, 1, 0));
		E(4, 2, 4, 0) = (one2p/2.0)*(E(3, 1, 4, 0) + 2*E(3, 2, 3, 0));
		E(4, 4, 3, 0) = (one2p/4.0)*(4*E(3, 3, 3, 0) + 3*E(3, 4, 2, 0));
		E(4, 3, 4, 0) = (one2p/4.0)*(3*E(3, 2, 4, 0) + 4*E(3, 3, 3, 0));
		E(4, 4, 4, 0) = one2p*(E(3, 3, 4, 0) + E(3, 4, 3, 0));
	}	
	
	return E;
}	

// Calculate the nuclear attraction between two primitives
// and a centre c, using the McMurchie Davidson scheme
double IntegralEngine::mmNucAttract(const PBF& u, const PBF& v, const Vector& ucoords, 
				  const Vector& vcoords, const Vector& ccoords) const
{
  double integral = 0.0; // To return the answer in

  // Get exponents, norms, and angular momenta
  int ulx = u.getLx(); int uly = u.getLy(); int ulz = u.getLz();
  int vlx = v.getLx(); int vly = v.getLy(); int vlz = v.getLz();
  double unorm = u.getNorm(); double uexp = u.getExponent();
  double vnorm = v.getNorm(); double vexp = v.getExponent();

  // Determine the maximum N needed for the auxiliary integrals
  // in each cartesian direction
  int Nx = ulx + vlx; int Ny = uly + vly; int Nz = ulz + vlz;
  int N = Nx + Ny + Nz;
  
  // Get the necessary values from getVals
  Vector vals;
  vals = getVals(uexp, vexp, ucoords, vcoords);
  
  // Calculate/retrieve the needed atomic separations
  double XPA = vals(2) - ucoords(0); double XPC = vals(2) - ccoords(0);
  double YPA = vals(3) - ucoords(1); double YPC = vals(3) - ccoords(1);
  double ZPA = vals(4) - ucoords(2); double ZPC = vals(4) - ccoords(2);
  double XPB = vals(2) - vcoords(0); double YPB = vals(3) - vcoords(1);
  double ZPB = vals(4) - vcoords(2); double p = vals(0);
  
  // Generate the expansion coefficients
  Tensor4 EX(N+1, N+1, N+1, 1);
  Tensor4 EY(N+1, N+1, N+1, 1);
  Tensor4 EZ(N+1, N+1, N+1, 1);
  
  EX = makeE(ulx+uly+ulz, vlx+vly+vlz, vals(8), p, XPA, XPB);
  EY = makeE(ulx+uly+ulz, vlx+vly+vlz, vals(9), p, YPA, YPB);
  EZ = makeE(ulx+uly+ulz, vlx+vly+vlz, vals(10), p, ZPA, ZPB);
  
  // Evaluate the Hermite integrals
  Tensor4 R(N+1, N+1, N+1, N+1, 0.0);

  double m2pn = -2.0*p; 
  double pRPC2 = p*(XPC*XPC + YPC*YPC + ZPC*ZPC);
  
  // Get the initial values from the Boys function
  Vector boysvals = boys(pRPC2, N, 0);
  for (int i = 0; i < N+1; i++)
  	R(i, 0, 0, 0) = std::pow(m2pn, i)*boysvals(i);
  	
  // Make the integrals, avoiding recursion
  if (N > 0){
  	R(N-1, 1, 0, 0) = XPC*R(N, 0, 0, 0);
  	R(N-1, 0, 1, 0) = YPC*R(N, 0, 0, 0);
  	R(N-1, 0, 0, 1) = ZPC*R(N, 0, 0, 0);
  }
  
  if (N > 1){
  	R(N-2, 1, 0, 0) = XPC*R(N-1, 0, 0, 0);
  	R(N-2, 0, 1, 0) = YPC*R(N-1, 0, 0, 0);
  	R(N-2, 0, 0, 1) = ZPC*R(N-1, 0, 0, 0);
  	R(N-2, 1, 1, 0) = YPC*R(N-1, 1, 0, 0);
  	R(N-2, 0, 1, 1) = ZPC*R(N-1, 0, 1, 0);
  	R(N-2, 1, 0, 1) = XPC*R(N-1, 0, 0, 1);
  	R(N-2, 2, 0, 0) = R(N-1, 0, 0, 0) + XPC*R(N-1, 1, 0, 0);
  	R(N-2, 0, 2, 0) = R(N-1, 0, 0, 0) + YPC*R(N-1, 0, 1, 0);
  	R(N-2, 0, 0, 2) = R(N-1, 0, 0, 0) + ZPC*R(N-1, 0, 0, 1);
  }
  
  if (N>2) {
  	R(N-3, 1, 0, 0) = XPC*R(N-2, 0, 0, 0);
  	R(N-3, 0, 1, 0) = YPC*R(N-2, 0, 0, 0);
  	R(N-3, 0, 0, 1) = ZPC*R(N-2, 0, 0, 0);
  	R(N-3, 1, 1, 0) = YPC*R(N-2, 1, 0, 0);
  	R(N-3, 0, 1, 1) = ZPC*R(N-2, 0, 1, 0);
  	R(N-3, 1, 0, 1) = XPC*R(N-2, 0, 0, 1);
  	R(N-3, 2, 0, 0) = R(N-2, 0, 0, 0) + XPC*R(N-2, 1, 0, 0);
  	R(N-3, 0, 2, 0) = R(N-2, 0, 0, 0) + YPC*R(N-2, 0, 1, 0);
  	R(N-3, 0, 0, 2) = R(N-2, 0, 0, 0) + ZPC*R(N-2, 0, 0, 1);
  	R(N-3, 1, 1, 1) = ZPC*R(N-2, 1, 1, 0);
  	R(N-3, 2, 1, 0) = R(N-2, 0, 1, 0) + XPC*R(N-2, 1, 1, 0);
  	R(N-3, 1, 2, 0) = R(N-2, 1, 0, 0) + YPC*R(N-2, 1, 1, 0);
  	R(N-3, 0, 2, 1) = R(N-2, 0, 0, 1) + YPC*R(N-2, 0, 1, 1);
  	R(N-3, 0, 1, 2) = R(N-2, 0, 1, 0) + ZPC*R(N-2, 0, 1, 1);
  	R(N-3, 2, 0, 1) = R(N-2, 0, 0, 1) + XPC*R(N-2, 1, 0, 1);
  	R(N-3, 1, 0, 2) = R(N-2, 1, 0, 0) + ZPC*R(N-2, 1, 0, 1);
  	R(N-3, 3, 0, 0) = 2*R(N-2, 1, 0, 0) + XPC*R(N-2, 2, 0, 0);
  	R(N-3, 0, 3, 0) = 2*R(N-2, 0, 1, 0) + YPC*R(N-2, 0, 2, 0);
  	R(N-3, 0, 0, 3) = 2*R(N-2, 0, 0, 1) + ZPC*R(N-2, 0, 0, 2);
  }
  
  if (N>3) {
  	R(N-4, 1, 0, 0) = XPC*R(N-3, 0, 0, 0);
  	R(N-4, 0, 1, 0) = YPC*R(N-3, 0, 0, 0);
  	R(N-4, 0, 0, 1) = ZPC*R(N-3, 0, 0, 0);
  	R(N-4, 1, 1, 0) = YPC*R(N-3, 1, 0, 0);
  	R(N-4, 0, 1, 1) = ZPC*R(N-3, 0, 1, 0);
  	R(N-4, 1, 0, 1) = XPC*R(N-3, 0, 0, 1);
  	R(N-4, 2, 0, 0) = R(N-3, 0, 0, 0) + XPC*R(N-3, 1, 0, 0);
  	R(N-4, 0, 2, 0) = R(N-3, 0, 0, 0) + YPC*R(N-3, 0, 1, 0);
  	R(N-4, 0, 0, 2) = R(N-3, 0, 0, 0) + ZPC*R(N-3, 0, 0, 1);
  	R(N-4, 1, 1, 1) = ZPC*R(N-3, 1, 1, 0);
  	R(N-4, 2, 1, 0) = R(N-3, 0, 1, 0) + XPC*R(N-3, 1, 1, 0);
  	R(N-4, 1, 2, 0) = R(N-3, 1, 0, 0) + YPC*R(N-3, 1, 1, 0);
  	R(N-4, 0, 2, 1) = R(N-3, 0, 0, 1) + YPC*R(N-3, 0, 1, 1);
  	R(N-4, 0, 1, 2) = R(N-3, 0, 1, 0) + ZPC*R(N-3, 0, 1, 1);
  	R(N-4, 2, 0, 1) = R(N-3, 0, 0, 1) + XPC*R(N-3, 1, 0, 1);
  	R(N-4, 1, 0, 2) = R(N-3, 1, 0, 0) + ZPC*R(N-3, 1, 0, 1);
  	R(N-4, 3, 0, 0) = 2*R(N-3, 1, 0, 0) + XPC*R(N-3, 2, 0, 0);
  	R(N-4, 0, 3, 0) = 2*R(N-3, 0, 1, 0) + YPC*R(N-3, 0, 2, 0);
  	R(N-4, 0, 0, 3) = 2*R(N-3, 0, 0, 1) + ZPC*R(N-3, 0, 0, 2);
  	R(N-4, 2, 1, 1) = R(N-3, 0, 1, 1) + XPC*R(N-3, 1, 1, 1);
  	R(N-4, 1, 2, 1) = R(N-3, 1, 0, 1) + YPC*R(N-3, 1, 1, 1);
  	R(N-4, 1, 1, 2) = R(N-3, 1, 1, 0) + ZPC*R(N-3, 1, 1, 1);
  	R(N-4, 2, 2, 0) = R(N-3, 0, 2, 0) + XPC*R(N-3, 1, 2, 0);
  	R(N-4, 0, 2, 2) = R(N-3, 0, 0, 2) + YPC*R(N-3, 0, 1, 2);
  	R(N-4, 2, 0, 2) = R(N-3, 2, 0, 0) + ZPC*R(N-3, 2, 0, 1);
  	R(N-4, 3, 1, 0) = 2*R(N-3, 1, 1, 0) + XPC*R(N-3, 2, 1, 0);
  	R(N-4, 3, 0, 1) = 2*R(N-3, 1, 0, 1) + XPC*R(N-3, 2, 0, 1);
  	R(N-4, 1, 3, 0) = 2*R(N-3, 1, 1, 0) + YPC*R(N-3, 1, 2, 0);
  	R(N-4, 1, 0, 3) = 2*R(N-3, 1, 0, 1) + ZPC*R(N-3, 1, 0, 2);
  	R(N-4, 0, 3, 1) = 2*R(N-3, 0, 1, 1) + YPC*R(N-3, 0, 2, 1);
  	R(N-4, 0, 1, 3) = 2*R(N-3, 0, 1, 1) + ZPC*R(N-3, 0, 1, 2);
  	R(N-4, 4, 0, 0) = 3*R(N-3, 2, 0, 0) + XPC*R(N-3, 3, 0, 0);
  	R(N-4, 0, 4, 0) = 3*R(N-3, 0, 2, 0) + YPC*R(N-3, 0, 3, 0);
  	R(N-4, 0, 0, 4) = 3*R(N-3, 0, 0, 2) + ZPC*R(N-3, 0, 0, 3);
  }	
  				
  // Sum these into the Cartesian integral
  for (int t = 0; t < N+1; t++){
  	for (int u = 0; u < N+1-t; u++){
  		for (int v = 0; v < N+1-t-u; v++){
  			integral += EX(t, ulx, vlx, 0)*EY(u, uly, vly, 0)*EZ(v, ulz, vlz, 0)*R(0, t, u, v);
  		}
  	}		
  }
  
  integral = ((2.0*M_PI)/p)*integral;
  return unorm*vnorm*integral;  
}				  


// Calculate the two-electron integrals over a shell
// quartet of basis functions, using the Obara-Saika
// horizontal recursion relations.
// Algorithm:
//      - Loop over primitive quartets
//         - Get the auxiliary integrals [u0 | w0]
//         - Contract these into the relevant cgbf integrals
//           of the form (m0|p0)
//      -Loop over cgbf quartets
//        - Use the horizontal recursion relations on the 
//          second electron to form (m0|pq)
//        - Sphericalise to (m0|cd)
//        - Use the horizontal recursion on first electron
//          to get (mn|cd)
//        - Sphericalise to (ab|cd)
Tensor4 IntegralEngine::twoe(Atom& A, Atom& B, Atom& C, Atom& D, 
			    int shellA, int shellB, int shellC, int shellD) const
{
  // Get the coordinates of each atom, number of prims in each shell,
  // no. of cgbfs in each shell 
  Vector cA; Vector cB; Vector cC; Vector cD;
  cA = A.getCoords(); cB = B.getCoords(); cC = C.getCoords(); cD = D.getCoords();
  int npA = A.getNShellPrims(shellA); int npB = B.getNShellPrims(shellB);
  int npC = C.getNShellPrims(shellC); int npD = D.getNShellPrims(shellD);
  Vector sA; Vector sB; Vector sC; Vector sD;
  sA = A.getShells(); sB = B.getShells(); sC = C.getShells(); sD = D.getShells();
  int ncA = sA(shellA); int ncB = sB(shellB); int ncC = sC(shellC); int ncD = sD(shellD);

  // Get the Lnums of the shells
  int LA = A.getShellBF(shellA, 0).getLnum();
  int LB = B.getShellBF(shellB, 0).getLnum();
  int LC = C.getShellBF(shellC, 0).getLnum();
  int LD = D.getShellBF(shellD, 0).getLnum();
  
  // Form a 4-tensor of 6-tensors to store prim ints and to contract into
  Ten4Ten6 prims(npA, npB, npC, npD);
  Ten4Ten6 contr(ncA, ncB, ncC, ncD);
  
  // Now we need to do all the calculations over primitive shell quartets
  PBF upbf; PBF vpbf; PBF wpbf; PBF xpbf;
  for (int u = 0; u < npA; u++){
    upbf = A.getShellPrim(shellA, u);

    for (int v = 0; v < npB; v++){
      vpbf = B.getShellPrim(shellB, v);

      for (int w = 0; w < npC; w++){
	wpbf = C.getShellPrim(shellC, w);

	for (int x = 0; x < npD; x++){
	  xpbf = D.getShellPrim(shellD, x);

	  // Calculate the primitive quartet integrals
	  prims.set(u, v, w, x, twoe(upbf, vpbf, wpbf, xpbf, cA, cB, cC, cD));
		} // End x-loop
      } // End w-loop
    } // End v-loop
  } // End u-loop


  // Contract prims into contr
  for (int a = 0; a < ncA; a++){
    Vector aplist; Vector acoeffs;
    aplist = A.getShellBF(shellA, a).getPrimList();
    acoeffs = A.getShellBF(shellA, a).getCoeffs();
    
    for (int b = 0; b < ncB; b++){
      Vector bplist; Vector bcoeffs;
      bplist = B.getShellBF(shellB, b).getPrimList();
      bcoeffs = B.getShellBF(shellB, b).getCoeffs();
      int blx = B.getShellBF(shellB, b).getLx();
      int bly = B.getShellBF(shellB, b).getLy();
      int blz = B.getShellBF(shellB, b).getLz();	
      
      for (int c = 0; c < ncC; c++){
	Vector cplist; Vector ccoeffs;
	cplist = C.getShellBF(shellC, c).getPrimList();
	ccoeffs = C.getShellBF(shellC, c).getCoeffs();
	
	for (int d = 0; d < ncD; d++){
	  Vector dplist; Vector dcoeffs;
	  dplist = D.getShellBF(shellD, d).getPrimList();
	  dcoeffs = D.getShellBF(shellD, d).getCoeffs();
	  int dlx = D.getShellBF(shellD, d).getLx();
	  int dly = D.getShellBF(shellD, d).getLy();
	  int dlz = D.getShellBF(shellD, d).getLz();
	  
	  contr(a, b, c, d).assign(blx+1, bly+1, blz+1, dlx+1, dly+1, dlz+1, 0.0);
	  for (int u = 0; u < aplist.size(); u++){
	    for (int v = 0; v < bplist.size(); v++){
	      for (int w = 0; w < cplist.size(); w++){
		for (int x = 0; x < dplist.size(); x++){
		  double cmult = acoeffs(u)*bcoeffs(v)*ccoeffs(w)*dcoeffs(x);
		  contr(a, b, c, d) = contr(a, b, c, d) + 
		    cmult*prims(aplist(u), bplist(v), cplist(w), dplist(x));
		}
	      }
	    }
	  }
	} // End of d-loop
      } // End of c-loop
    } // End of b-loop
  } // End of a-loop									
  
  // Done with prims
  prims.clean();	
  
  // Get atomic separations
  double XAB, YAB, ZAB, XCD, YCD, ZCD;
  XAB = cA(0)-cB(0); YAB = cA(1)-cB(1); ZAB = cA(2)-cB(2);
  XCD = cC(0)-cD(0); YCD = cC(1)-cD(1); ZCD = cC(2)-cD(2);
  
  // We now have contracted integrals of the form (m0|p0) sitting in 
  // the contr ten4ten6. First we transform these to (m0|pq) by the
  // horizontal recursion relation.
  int nlx, nly, nlz, qlx, qly, qlz;
  for (int m = 0; m < ncA; m++){
    for (int n = 0; n < ncB; n++){
      // Get angular momenta
      nlx = B.getShellBF(shellB, n).getLx();
      nly = B.getShellBF(shellB, n).getLy();
      nlz = B.getShellBF(shellB, n).getLz();
      
      for (int p = 0; p < ncC; p++){
	for (int q = 0; q < ncD; q++){
	  // Get angular momenta
	  qlx = D.getShellBF(shellD, q).getLx();
	  qly = D.getShellBF(shellD, q).getLy();
	  qlz = D.getShellBF(shellD, q).getLz();
	  // Increment q
  	  for (int vx = 0; vx < nlx+1; vx++){
	    for (int vy = 0; vy < nly+1; vy++){
	      for (int vz = 0; vz < nlz+1; vz++){
		for (int inc = 1; inc < qlx+1; inc++){		
		  // Start with the x-coordinate
		  for (int xy = 0; xy < qly+1; xy++){
		    for (int xz = 0; xz < qlz+1; xz++){
		      for (int xx = 0; xx < qlx-inc+1; xx++){
			contr(m, n, p, q)(vx, vy, vz, xx, xy, xz) = 
			  contr(m, n, p, q)(vx, vy, vz, xx+1, xy, xz) 
			  + XCD*contr(m, n, p, q)(vx, vy, vz, xx, xy, xz);
		      }
		    }
		  }
		}
	  	
		// Then the y-coordinate
		for (int inc = 1; inc < qly+1; inc++){
		  for (int xz = 0; xz < qlz+1; xz++){
		    for (int xy = 0; xy < qly-inc+1; xy++){
		      contr(m, n, p, q)(vx, vy, vz, 0, xy, xz) = 
			contr(m, n, p, q)(vx, vy, vz, 0, xy+1, xz) 
			+ YCD*contr(m, n, p, q)(vx, vy, vz, 0, xy, xz);
		    }
		  }
		}
		
		// And finally the z-coordinate	 
		for (int inc = 1; inc < qlz+1; inc++){
		  for (int xz = 0; xz < qlz-inc+1; xz++){
		    contr(m, n, p, q)(vx, vy, vz, 0, 0, xz) = 
		      contr(m,n,p,q)(vx,vy,vz,0,0, xz+1) + ZCD*contr(m,n,p,q)(vx,vy,vz,0,0,xz);
		  }
		}
	      } // End of vz-loop
	    } // End of vy-loop					
	  } // End of vx-loop 
	} // End of q-loop
      } // End of p-loop
    } // End of n-loop
  } // End of m-loop

  // The integrals are all now of the form (m0|pq), and the second electron is ready to be
  // transformed to the spherical harmonic basis. 
  
  // Get the number of transformed basis functions on C/D
  int spherC; int spherD;
  switch(LC){
  case 2: { spherC = 5*(sC(shellC)/6); break; }
  case 3: { spherC = 7*(sC(shellC)/10); break; }
  default: spherC = ncC;
  }
  switch(LD){
  case 2: { spherD = 5*(sD(shellD)/6); break; }
  case 3: { spherD = 7*(sD(shellD)/10); break; }
  default: spherD = ncD;
  }
  
  // Make a list of mnums for each
  Vector cmnums(2*LC+1); Vector dmnums(2*LD+1);
  
  // Maps 0 -> 0, 1 -> -1, 2 -> 1, 3 -> -2, 4 -> 2, and so on.
  for (int cm = 0; cm < cmnums.size(); cm++)
    cmnums[cm] = (1-2*(cm%2))*((cm+1)/2);

  if (LC == 2) { 
    cmnums[1] = -2; cmnums[3] = 2; cmnums[4] = -1;
  } else if (LC == 1){ cmnums[0] = 1; cmnums[2] = 0; }

  for (int dm = 0; dm < dmnums.size(); dm++)
  	dmnums[dm] = (1-2*(dm%2))*((dm+1)/2);	
  
  if (LD == 2){
    dmnums[1] = -2; dmnums[3] = 2; dmnums[4] = -1;
  } else if (LD== 1){ dmnums[0] = 1; dmnums[2] = 0; }

  // Declare appropriately sized C and D transformation matrix
  Matrix tMat1(spherC, ncC, 0.0); Matrix tMat2(spherD, ncD, 0.0);
  
  // Build the C trans matrix
  int jinc; 
  switch(LC){
  	case 1: { jinc = 3; break; }
  	case 2: { jinc = 6; break; }
  	case 3: { jinc = 10; break; }
  	default: jinc = 1;
  }
  int j = 0; // Column counter
  int mod = 2*LC+1;
  int mcount = 0;
  for(int i = 0; i < spherC; i++){  	
    // Get the coefficients
    formTransMat(tMat1, i, j, LC, cmnums(i%mod));
    if (mcount == 2*LC){ // Increment j by a suitable amount
      j+=jinc; mcount = 0;
    } else { mcount++; }
  }
  
  // Build the D trans matrix
  mcount = 0;
  switch(LD){
  case 1: { jinc = 3; break; }
  case 2: { jinc = 6; break; }
  case 3: { jinc = 10; break; }
  default: jinc = 1;
  }
  j = 0; mod = 2*LD+1;

  for(int i = 0; i < spherD; i++){  	
    // Get the coefficients
    formTransMat(tMat2, i, j, LD, dmnums(i%mod));
    if (mcount == 2*LD){ // Increment j by a suitable amount
      j+=jinc; mcount = 0;
    } else { mcount++; }
  }
  
  // Transform contr into halfspher
  Ten4Ten4 halfspher(ncA, ncB, spherC, spherD);
  for (int m = 0; m < ncA; m++){
    for (int n = 0; n < ncB; n++){
      nlx = B.getShellBF(shellB, n).getLx();
      nly = B.getShellBF(shellB, n).getLy();
      nlz = B.getShellBF(shellB, n).getLz();
      
      // Resize the relevant halfspher tensors
      for (int c = 0; c < spherC; c++){
	for (int d = 0; d < spherD; d++){
	  halfspher(m, n, c, d).assign(nlx+1, nly+1, nlz+1, 1, 0.0);
	}
      }
      // Construct and transform all the pmats
      for (int x = 0; x < nlx+1; x++){
	for (int y = 0; y < nly+1; y++){
	  for (int z = 0; z < nlz+1; z++){
	    Matrix temp(ncC, ncD, 0.0);
	    for (int p = 0; p < ncC; p++){
	      for (int q = 0; q < ncD; q++) {
		temp(p, q) = contr(m, n, p, q)(x, y, z, 0, 0, 0);
	      }
	    }
	    
	    // Transform and copy into halfsher
	    for (int c = 0; c < spherC; c++){
	      for (int d = 0; d < spherD; d++){
		for (int p = 0; p < ncC; p++){
		  for (int q = 0; q < ncD; q++){
		    halfspher(m, n, c, d)(x, y, z, 0) += tMat1(c, p)*temp(p, q)*tMat2(d, q);
		  }
		}		
	      }
	    } // End of copy
	  } // End of z-loop
	} // End of y-loop
      } // End of x-loop
      
    } // End of n-loop
  } // End of m-loop
  
  // Get rid of contr, as everything is now in halfspher
  contr.clean();
  
  // We now have all integrals of the form (m0|cd).
  // Move on to the second horizontal recursion step.
  for (int m = 0; m < ncA; m++){
    for (int n = 0; n < ncB; n++){
      nlx = B.getShellBF(shellB, n).getLx();
      nly = B.getShellBF(shellB, n).getLy();
      nlz = B.getShellBF(shellB, n).getLz();
      
      for (int c = 0; c < spherC; c++){
	for (int d = 0; d < spherD; d++){
	  
	  // Increment in the x-index
	  for (int inc = 1; inc < nlx+1; inc++){
	    for (int ny = 0; ny < nly+1; ny++){
	      for (int nz = 0; nz < nlz+1; nz++){
		for (int nx = 0; nx < nlx-inc+1; nx++){
		  halfspher(m,n,c,d)(nx,ny,nz,0) = halfspher(m,n,c,d)(nx+1,ny,nz,0)
		    +XAB*halfspher(m,n,c,d)(nx,ny,nz,0);
		}
	      }
	    }
	  }					
	  
	  // Increment in the y-index
	  for (int inc = 1; inc < nly+1; inc++){
	    for (int nz = 0; nz < nlz+1; nz++){
	      for (int ny = 0; ny < nly-inc+1; ny++){
		halfspher(m,n,c,d)(0, ny, nz, 0) = halfspher(m,n,c,d)(0,ny+1,nz,0)
		  +YAB*halfspher(m,n,c,d)(0,ny,nz,0);
	      }
	    }			
	  }
	  
	  // Finally, increment in the z-index
	  for (int inc = 1; inc < nlz+1; inc++){
	    for (int nz = 0; nz < nlz-inc+1; nz++){
	      halfspher(m,n,c,d)(0, 0, nz, 0) = halfspher(m,n,c,d)(0,0,nz+1,0)
		+ ZAB*halfspher(m,n,c,d)(0,0,nz,0);
	    }		
	  }
	} // End of d-loop
      } // End of c-loop
    } // End of n-loop
  } // End of m-loop
  
  // All integrals are now (mn|cd), stored in the first element of each tensor
  // Only thing left to	do is to sphericalise the first electron
  
  // Get the number of transformed basis functions on A/B
  int spherA; int spherB;
  switch(LA){
  case 2: { spherA = 5*(sA(shellA)/6); break; }
  case 3: { spherA = 7*(sA(shellA)/10); break; }
  default: spherA = ncA;
  }
  switch(LB){
  case 2: { spherB = 5*(sB(shellB)/6); break; }
  case 3: { spherB = 7*(sB(shellB)/10); break; }
  default: spherB = ncB;
  }
  
  // Make a list of mnums for each
  cmnums.resize(2*LA+1); dmnums.resize(2*LB+1);
  
  // Maps 0 -> 0, 1 -> -1, 2 -> 1, 3 -> -2, 4 -> 2, and so on.
  for (int cm = 0; cm < cmnums.size(); cm++)
    cmnums[cm] = (1-2*(cm%2))*((cm+1)/2);
  
  if (LA == 2) { cmnums[1] = -2; cmnums[3] = 2; cmnums[4] = -1; }
  else if (LA== 1){ cmnums[0] = 1; cmnums[2] = 0; }

  for (int dm = 0; dm < dmnums.size(); dm++)
    dmnums[dm] = (1-2*(dm%2))*((dm+1)/2);	
  
  if (LB == 2) { dmnums[1] = -2; dmnums[3] = 2; dmnums[4] = -1; }
  else if (LB== 1){ dmnums[0] = 1; dmnums[2] = 0; }

  // Resize the tMats
  tMat1.assign(spherA, ncA, 0.0); tMat2.assign(spherB, ncB, 0.0);
  
  // Build the A trans matrix
  switch(LA){
  case 1: { jinc = 3; break; }
  case 2: { jinc = 6; break; }
  case 3: { jinc = 10; break; }
  default: jinc = 1;
  }
  j = 0; // Column counter
  mod = 2*LA+1;
  mcount = 0;
  for(int i = 0; i < spherA; i++){  	
    // Get the coefficients
    formTransMat(tMat1, i, j, LA, cmnums(i%mod));
    if (mcount == 2*LA){ // Increment j by a suitable amount
      j+=jinc; mcount = 0;
    } else { mcount++; }
  }
  
  // Build the B trans matrix
  mcount = 0;
  switch(LB){
  case 1: { jinc = 3; break; }
  case 2: { jinc = 6; break; }
  case 3: { jinc = 10; break; }
  default: jinc = 1;
  }
  j = 0; mod = 2*LB+1;
  for(int i = 0; i < spherB; i++){  	
    // Get the coefficients
    formTransMat(tMat2, i, j, LB, dmnums(i%mod));
    if (mcount == 2*LB){ // Increment j by a suitable amount
      j+=jinc; mcount = 0;
    } else { mcount++; }
  }
  
  // Transform the integrals into the return vector retInts
  Tensor4 retInts(spherA, spherB, spherC, spherD, 0.0);
  for (int c = 0; c < spherC; c++){	
    for (int d = 0; d < spherD; d++){
      
      Matrix temp(ncA, ncB, 0.0);
      // Construct temp and transform
      for (int m = 0; m < ncA; m++){
	for (int n = 0; n < ncB; n++){
	  temp(m, n) = halfspher(m, n, c, d)(0, 0, 0, 0);
	}
      }		 
      
      // Copy into retInts
      for (int a = 0; a < spherA; a++){
	for (int b = 0; b < spherB; b++){
	  for (int m = 0; m < ncA; m++){
	    for (int n = 0; n < ncB; n++){
	      retInts(a,b,c,d) += tMat1(a, m)*temp(m, n)*tMat2(b, n);
	    }
	  }		
	}
      } // End of retInts copy
      
    } // End of d-loop
  } // End of c-loop  
  
  // Get rid of halfspher
  halfspher.clean();
  
  // Integrals are now all (ab|cd) and are stored in retInts
  return retInts;
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
//   Return these as a Tensor ready for contraction, 
//   sphericalisation, and then horizontal recurrence.
Tensor6 IntegralEngine::twoe(const PBF& u, const PBF& v, const PBF& w,
			    const PBF& x, const Vector& ucoords,
			    const Vector& vcoords, const Vector& wcoords,
			    const Vector& xcoords) const
{
  // Extract all the necessary data from the PBFs
  Vector pvals; Vector qvals;
  pvals = getVals(u.getExponent(), v.getExponent(), ucoords, vcoords);
  qvals = getVals(w.getExponent(), x.getExponent(), wcoords, xcoords);
  
  // Unpack, and calculate distances, exponents, and multipliers.
  double p = pvals(0); double q = qvals(0); double alpha = (p*q)/(p+q);
  double XPA = pvals(2) - ucoords(0); double XPQ = pvals(2) - qvals(2);
  double YPA = pvals(3) - ucoords(1); double YPQ = pvals(3) - qvals(3);
  double ZPA = pvals(4) - ucoords(2); double ZPQ = pvals(4) - qvals(4);
  double XAB = pvals(5); double YAB = pvals(6); double ZAB = pvals(7);
  double XCD = qvals(5); double YCD = qvals(6); double ZCD = qvals(7);
  double K = pvals(8)*pvals(9)*pvals(10)*qvals(8)*qvals(9)*qvals(10);
  double zeromult = 2.0*K*M_PI*M_PI*std::sqrt(M_PI/(p+q))/(p*q);
  double RPQ2 = XPQ*XPQ + YPQ*YPQ + ZPQ*ZPQ;
  double ap = alpha/p; double one2p = 1.0/(2.0*p); double one2q = 1.0/(2.0*q);
  double poq = p/q; double norms = u.getNorm()*v.getNorm()*w.getNorm()*x.getNorm();
  zeromult *= norms;
  
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
  Vector boysvals(L+1, 0.0); 
  
  // Make a tensor to store results in
  Tensor4 aux(L+1, Nx+1, Ny+1, Nz+1, 0.0);

  // First calculate all boys function values
  boysvals = boys(alpha*RPQ2, L, 0);
  
  // Then convert 
  for (int i = 0; i < L+1; i++)
    aux(i, 0, 0, 0) = zeromult*boysvals(i);

  // Now we need to use the vertical recurrence relation to increment
  // the first index, one cartesian direction at a time
  if ( Nx > 0 ){
  for (int n = L-1; n > -1; n--){
    // Calculate first increment
    aux(n, 1, 0, 0) = XPA*aux(n, 0, 0, 0) - ap*XPQ*aux(n+1, 0, 0, 0);
    
    // Now increment as high as needed
    int kmax = (L-n+1 > Nx+1 ? Nx+1 : L-n+1);
    for (int k = 2; k < kmax; k++){
      aux(n, k, 0, 0) = XPA*aux(n, k-1, 0, 0) - ap*XPQ*aux(n+1, k-1, 0, 0) + (k-1)*one2p*(aux(n, k-2, 0, 0) - ap*aux(n+1, k-2, 0, 0));
    }
  }
  }

  // for the Nx+1 columns truncated at n = L-Nx, increment
  // the y-direction on u
  
  for (int m = 0; m < Nx+1; m++){

    // Now increment the y-direction up to Ny
    if (Ny > 0) {
    for (int n = L-Nx-1; n > -1; n--){
      // Calculate the first increment
      aux(n, m, 1, 0) = YPA*aux(n, m, 0, 0) - ap*YPQ*aux(n+1, m, 0, 0);

      int kmax = (L-Nx-n+1 > Ny+1 ? Ny+1 : L-Nx-n+1);
      for (int k = 2; k < kmax; k++){
		aux(n, m, k, 0) = YPA*aux(n, m, k-1, 0) - ap*YPQ*aux(n+1, m, k-1, 0)
		+ (k-1)*one2p*(aux(n, m, k-2, 0) - ap*aux(n+1, m, k-2, 0));
      }
    }
    }
    
  }

  // Next increment in the z-direction
  for (int m = 0; m < Nx+1; m++){
  	for (int n = 0; n < Ny+1; n++){

	    // Now increment
    	for (int row = Nz-1; row > -1; row--){
      		// Calculate first increment
     		 aux(row, m, n, 1) = ZPA*aux(row, m, n, 0) - ap*ZPQ*aux(row+1, m, n, 0);
      
     		 for (int k = 2; k < Nz-row+1; k++){
				aux(row, m, n, k) = ZPA*aux(row, m, n, k-1) - ap*ZPQ*aux(row+1, m, n, k-1) 
				+ (k-1)*one2p*(aux(row, m, n, k-2) - ap*aux(row+1, m, n, k-2));
     		 }
    	}  
    }	  
  }

  // We now have all the integrals [u0|00] stored in aux
  // for u in [0, L]. 
  
  // Next step is then to use the electron-transfer recurrence relation
  // to transfer cartesian powers from u to w (i.e. from elec 1 to elec 2)
  // We need [u0|w0] for u in [Lu, Lu+Lv] and w in [Lw, Lw+Lx]
	
  // Some premultipliers	
  double multX = -(v.getExponent()*XAB + x.getExponent()*XCD)/q;
  double multY = -(v.getExponent()*YAB + x.getExponent()*YCD)/q;
  double multZ = -(v.getExponent()*ZAB + x.getExponent()*ZCD)/q;
  
  // Get the components of the angular momenta
  int wlx = w.getLx(); int xlx = x.getLx(); int ulx = u.getLx(); int vlx = v.getLx();
  int wly = w.getLy(); int xly = x.getLy(); int uly = u.getLy(); int vly = v.getLy();
  int wlz = w.getLz(); int xlz = x.getLz(); int ulz = u.getLz(); int vlz = v.getLz();

  // Make a tensor for the calculations
  Tensor6 newAux(Nx+1, Ny+1, Nz+1, wlx+xlx+1, wly+xly+1, wlz+xlz+1, 0.0);
  
  // Copy aux into newAux
  for (int m = 0; m < Nx+1; m++){
    for (int n = 0; n < Ny+1; n++){
      for (int r = 0; r < Nz+1; r++){
	newAux(m, n, r, 0, 0, 0) = aux(0, m, n, r);
      }
    }		
  }
  
  // Get rid of aux
  aux.resize(0,0,0,0);
  
  // Increment the x-index
  if(wlx+xlx>0){
    // Do the first increment of first bit
    for (int n = 0; n < Ny+1; n++){
      for (int r = 0; r < Nz+1; r++){
	newAux(0, n, r, 1, 0, 0) = multX*newAux(0, n, r, 0, 0, 0) -poq*newAux(1, n, r, 0, 0, 0);
      }
    }
    
    // And the rest of first increment
    for (int m = 1; m < Nx; m++){
      for (int n = 0; n < Ny+1; n++){
	for (int r = 0; r < Nz+1; r++){
	  newAux(m, n, r, 1, 0, 0) = multX*newAux(m, n, r, 0, 0, 0) + m*one2q*newAux(m-1, n, r, 0, 0, 0)
	    -poq*newAux(m+1, n, r, 0, 0, 0);
	}	
      }	
    }
    
    // And the rest of the increments
    for (int inc = 2; inc < wlx+xlx+1; inc++){
      // Do the first bit
      for (int n = 0; n < Ny+1; n++){
	for (int r = 0; r < Nz+1; r++){
	  newAux(0, n, r, inc, 0, 0) = multX*newAux(0, n, r, inc-1, 0, 0) 
	    +(inc-1)*one2q*newAux(0, n, r, inc-2, 0, 0) - poq*newAux(1, n, r, inc-1, 0, 0);
	}
      }
      
      // And the rest
      for (int m = 1; m < Nx-inc+1; m++){
	for (int n = 0; n < Ny+1; n++){
	  for (int r = 0; r < Nz+1; r++){
	    newAux(m, n, r, inc, 0, 0) = multX*newAux(m, n, r, inc-1, 0, 0)
	      +m*one2q*newAux(m-1, n, r, inc-1, 0, 0) + (inc-1)*one2q*newAux(m, n, r, inc-2, 0, 0)
	      -poq*newAux(m+1, n, r, inc-1, 0, 0);
	  }
	}
      }
    } // End of further increments
  } // End of x-if-loop
  
  // Increment the y-index
  if (wly+xly>0){
    for (int m = ulx; m < ulx+vlx+1; m++){
      for (int s = wlx; s < wlx+xlx+1; s++){
	// Do the first increment of the first bit
	for (int r = 0; r < Nz+1; r++){
	  newAux(m, 0, r, s, 1, 0) = multY*newAux(m, 0, r, s, 0, 0) -poq*newAux(m, 1, r, s, 0, 0);
	}
  	
	// And the rest of the first increment
	for (int n = 1; n < Ny; n++){
	  for (int r = 0; r < Nz+1; r++){
	    newAux(m, n, r, s, 1, 0) = multY*newAux(m, n, r, s, 0, 0) + n*one2q*newAux(m, n-1, r, s, 0, 0)
	      -poq*newAux(m, n+1, r, s, 0, 0);
	  }
	}
  	
	// And the rest of the increments
	for (int inc = 2; inc < wly+xly+1; inc++){
	  // Do the first bit
	  for (int r = 0; r < Nz+1; r++){
	    newAux(m, 0, r, s, inc, 0) = multY*newAux(m, 0, r, s, inc-1, 0)
	      + (inc-1)*one2q*newAux(m, 0, r, s, inc-2, 0) -poq*newAux(m, 1, r, s, inc-1, 0);
	  }
	  
	  // And the rest
	  for (int n = 1; n < Ny-inc+1; n++){
	    for (int r = 0; r < Nz+1; r++){
	      newAux(m, n, r, s, inc, 0) = multY*newAux(m, n, r, s, inc-1, 0)
		+n*one2q*newAux(m, n-1, r, s, inc-1, 0) + (inc-1)*one2q*newAux(m, n, r, s, inc-2, 0)
		-poq*newAux(m, n+1, r, s, inc-1, 0);
	    }
	  }
	} // End of increment loop				
      } // End of s-loop	 			
    } // End of m-loop
  } // End of y-if-loop
  
  // Finally, increment the z-index
  if(wlz+xlz>0){
    for (int m = ulx; m < ulx+vlx+1; m++){
      for (int n = uly; n < uly+vly+1; n++){
	for (int s = wlx; s < wlx+xlx+1; s++){
	  for (int t = wly; t < wly+xly+1; t++){
	    // Do the first increment of first element
	    newAux(m, n, 0, s, t, 1) = multZ*newAux(m, n, 0, s, t, 0) -poq*newAux(m, n, 1, s, t, 0);
	    
	    // Do the rest of first increment
	    for (int r = 1; r < Nz; r++){
	      newAux(m, n, r, s, t, 1) = multZ*newAux(m, n, r, s, t, 0) 
		+r*one2q*newAux(m, n, r-1, s, t, 0) - poq*newAux(m, n, r+1, s, t, 0);
	    }
	    
	    // And the rest of the increments
	    for (int inc = 2; inc < wlz+xlz+1; inc++){
	      // Do the first element
	      newAux(m, n, 0, s, t, inc) = multZ*newAux(m, n, 0, s, t, inc-1)
		+(inc-1)*one2q*newAux(m, n, 0, s, t, inc-2) - poq*newAux(m, n, 1, s, t, inc-1);
	      
	      // And the rest
	      for (int r = 1; r < Nz-inc+1; r++){
		newAux(m, n, r, s, t, inc) = multZ*newAux(m, n, r, s, t, inc-1)
		  +r*one2q*newAux(m, n, r-1, s, t, inc-1) + (inc-1)*one2q*newAux(m, n, r, s, t, inc-2)
		  -poq*newAux(m, n, r+1, s, t, inc-1);
	      }			
	    } // End of increment loop
	  } // End of t-loop
	} //End of s-loop
      } // End of n-loop
    } // End of m-loop
  } // End of z-if-loop
  
  // All necessary integrals of the form [u0|w0] have now been calculated
  // and are stored in newAux. Transfer them to a more suitably sized tensor
  // for returning.
  
  Tensor6 retInts(vlx+1, vly+1, vlz+1, xlx+1, xly+1, xlz+1, 0.0);
  for (int m = ulx; m < ulx+vlx+1; m++){
    for (int n = uly; n < uly+vly+1; n++){
      for (int r = ulz; r < ulz+vlz+1; r++){
	for (int s = wlx; s < wlx+xlx+1; s++){
	  for (int t = wly; t < wly+xly+1; t++){
	    for (int z = wlz; z < wlz+xlz+1; z++){
	      retInts(m-ulx, n-uly, r-ulz, s-wlx, t-wly, z-wlz) = newAux(m, n, r, s, t, z);
	    }	
	  }
	}
      }
    }
  }						 						
  
  return retInts;
}

