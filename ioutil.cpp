/*
 *
 *   PURPOSE: Implements ioutil.hpp, a collection of i/o utility classes and functions.
 *
 *   DATE           AUTHOR           CHANGES
 *   ==================================================================
 *   30/08/15       Robert Shaw      Original code.
 *
 */

// Includes
#include "vector.hpp"
#include "basis.hpp"
#include "error.hpp"
#include <algorithm>

// Class FileReader implementation

// Destructor
FileReader::~FileReader
{
  if ( natoms!= 0 ) {
    delete[] geometry;
  }
}

// Match a token to a particular command
int FileReader::findToken(std::string t)
{
  // Get rid of extraneous spaces
  t.erase(std::remove(t.begin(), t.end(), ' '), t.end());
  
  // make t lowercase
  std::transform(t.begin(), t.end(), t.begin(), ::tolower);

  // Match it to a command
  int rval = 0;
  if (t == "charge") { rval = 1; }
  else if (t == "multiplicity") { rval = 2; }
  else if (t == "maxiter") { rval = 3; }
  else if (t == "precision") { rval = 4; }
  else if (t == "basis") { rval = 5; }
  else if (t == "geom") { rval = 6; }
  else if (t == "geomend") { rval = 7; }
  
  return rval;
}

// Read in the parameters from the input file
void FileReader::readParameters()
{
  // Set default values, in case none specified
  charge = 0;
  multiplicity = 1;
  maxiter = 50;
  precision = 1e-12;
  geomstart = 0; geomend = 0;

  // Read line by line and parse
  std::string line, token;
  int pos;
  int linecount = 0;
  while (std::getline(input, line)) {
    // Erase any comments
    pos = line.find('!');
    line.erase(pos, line.end());
    
    // Tokenise
    pos = line.find(',');
    if(pos != std::string::npos){
      token = line.substr(line.begin(), pos);
      // Match token to command
      switch(findToken(token)){
      case 1: { // Charge
	charge = std::stoi(line.substr(pos, line.end()));
	break;
      }
      case 2: { // Multiplicity
	multiplicity = std::stoi(line.substr(pos, line.end()));
	break;
      }
      case 3: { // MaxIter
	maxiter = std::stoi(line.substr(pos, line.end()));
	break;
      }
      case 4: { // Precision
	precision = std::stod(line.substr(pos, line.end()));
	break;
      }
      case 5: { // Basis
	line.erase(line.begin(), pos+1);
	// Get rid of excess spaces
	line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
	basis = line;
	break;
      }
      case 6: { // GeomStart
	geomstart = linecount + 1;
	break;
      }
      case 7: { // GeomEnd
	geomend = linecount;
	break;
      }
      default: { // Unkown command issued
	throw(Error("READIN", "Command " + token + " not found."));
	break;
      }
      }
    }
    linecount++;
  }
  natoms = geomend - geomstart;
}

void FileReader::readGeometry()
{
  // The array size will be geomend-geomstart
  if (natoms != 0){
    geometry = new std::string[natoms];
   
    // Rewind to beginning of file
    input.clear();
    input.seekg(0, std::ios::beg);
    
    // Copy the geometry into the array
    for (int l = 0; l < geomend; l++){
      std::getline(input, line);
      if (l > geomstart - 1){
	geometry[l - geomstart] = line;
      }
    }
  } else {
    throw(Error("READGM", "Geometry not found.")); 
  }
}

// Implement class BasisReader

// Open the basis file that contains the correct basis functions
void BasisReader::openFile(int q)
{
  // Find which row of elements q is in
  // Note that the first row is taken here to be H-Ne
  std::string row = "first";
  if (q > 10 && q < 19){ row = "second"; } // Na - Ar
  else if (q > 18 && q < 37) { row = "third"; } // K - Kr
  else if (q > 36 && q < 55) { row = "fourth"; } // Rb - Xe
  else if (q > 54 && q < 87) { row = "fifth"; } // Cs - Rn
  else if (q > 86) { row = "sixth"; }
  std::string filename = "basissets/";
  filename += name; filename += row; filename += ".basis";
  // Open file, read only
  input.open(filename, std::ifstream::in);
}

void BasisReader::closeFile()
{
  if(input.is_open()){
    input.close();
  }
}

// Read in the number of contracted gaussian basis functions
// associated with an atom of atomic number q in the basis set
int BasisReader::readNbfs(int q)
{
  int nbfs = 0;
  Vector shells = readShells(q);
  for (int i = 0; i < shells.size(); i++){
  	nbfs += shells(i);
  }
  return nbfs;
}

// Read in the ith basis function for atom q from the basis file
BF BasisReader::readBF(int q, int i)
{
  // Open the file
  openFile(q);

  std::string delim = ",";
  int l1 = 0, l2 = 0, l3 = 0; // Angular momenta
  Vector c; Vector e; // Coeffs and exps

  // Check it's open
  if (input.is_open()){
    // Parse
    std::string aname = getAtomName(q);
    aname += " "; aname += delim;
    std::size_t position;
    int bfcount = 0;
    std::string line, shell;
    std::getline(input, line);
    
    // Main loop
    while(!input.eof() && bfcount != i){
      position = line.find(aname);
      if (position != std::string::npos){ // Atom type found
	// Get the shell type
	position = line.find(delim);
	shell = line.substr(0, position);
	
	// Calculate total no. of exponents
	int nexps = std::count(line.begin(), line.end(), ',') - 1;
	// Copy in
	Vector tempexps(nexps);
	// Get rid of shell declaration
	line.erase(0, position+delim.length());
	position = line.find(delim);
	std::string temp2;
	int counter = 0;
	
	// Now front bits have been removed, copy in exps 1 by 1
	while (position != std::string::npos) {
	  line.erase(0, position+delim.length());
	  position = line.find(delim);
	  temp2 = line.substr(0, position); // Get the exponent
	  // Get rid of extraneous spaces
	  temp2.erase(std::remove(temp2.begin(), temp2.end(), ' '), temp2.end()); 
	  tempexps[counter] = std::stod(temp2);
	  counter++;
	}
	
	// How many functions does this shell have?
	// N.b. cartesian not spherical gaussians
	int lmult = 1;
	if(shell == "p") { lmult = 3; }
	else if (shell == "sp") { lmult = 4; }
	else if (shell == "d") { lmult = 6; }
	else if (shell == "f") { lmult = 10; }
	
	// Iterate through bfs to find right one
	int sublmult;
	std::getline(input, line);
	while (bfcount!=i && line.at(0) == 'c'){
		sublmult = 1;
	  while (bfcount != i && sublmult < lmult+1){
	    bfcount++;
	    sublmult++;
	  }
	  if (bfcount!=i){
	    std::getline(input, line);
	  }
	}
	
	if (bfcount == i) { // Found it 
	  // Find out which exponents to use
	  position = line.find(delim);
	  
	  // Get rid of exponent number declaration
	  line.erase(0, position+delim.length());
	  position = line.find(delim);
	  temp2 = line.substr(0, position);
	  int p = temp2.find('.');
	  int start = std::stoi(temp2.substr(0, p));
	  int end = std::stoi(temp2.substr(p+1, temp2.length()));
	  
	  // Now resize e and copy in
	  e.resize(end - start + 1);
	  for (int j = 0; j < end - start + 1; j++){
	    e[j] = tempexps(start+j-1); // Take account of zero-indexing
	  }
	  
	  // Get the coeffs
	  c.resize(end - start + 1);
	  for (int j = 0; j < end - start + 1; j++){
	    line.erase(0, position+delim.length());
	    position = line.find(delim);
	    temp2 = line.substr(0, position);
	    temp2.erase(std::remove(temp2.begin(), temp2.end(), ' '), temp2.end());
	    c[j] = std::stod(temp2);
	  }
	  
	  // Work out the lnums
	  // sublmult is 1 too high
	  sublmult--;
	  switch(lmult){
	  case 1: { // s type
	    l1 = l2 = l3 = 0;
	    break;
	  }
	  case 3: { // p type
	    switch(sublmult){
	    case 1: { //px 
	      l1 = 1; l2 = l3 = 0;
	      break;
	    }
	    case 2: { //py
	      l1 = l3 = 0; l2 = 1;
	      break;
	    }
	    case 3: { //pz
	      l1 = l2 = 0; l3 = 1;
	      break;
	    }
	    }
	    break;
	  }
	  case 4: { // sp type
	    switch(sublmult){
	    case 1:{ // s
	      l1 = l2 = l3 = 0;
	      break;
	    }
	    case 2:{ // px
	      l1 = 1; l2 = l3 = 0;
	      break;
	    }
	    case 3:{ // py
	      l1 = l3 = 0; l2 = 1;
	      break;
	    }
	    case 4:{ //pz
	      l1 = l2 = 0; l3 = 1;
	      break;
	    }
	    }
	    break;
	  }
	  case 6: { // d type
	    switch(sublmult){
	    case 1:{ // dzz
	      l3 = 2; l1 = l2 = 0;
	      break;
	    }
	    case 2:{ // dxy
	      l1 = l2 = 1; l3 = 0;
	      break;
	    }
	    case 3:{ // dxz
	      l1 = l3 = 1; l2 = 0;
	      break;
	    }
	    case 4:{ // dyz
	      l1 = 0; l2 = l3 = 1;
	      break;
	    }
	    case 5:{ // dxx 
	      l1 = 2; l2 = l3 = 0;
	      break;
	    }
	    case 6:{ // dyy
	      l2 = 2; l1 = l3 = 0;
	      break;
	    }
	    }
	    break;
	  }
	  case 10:{ // f type
	    switch(sublmult){
	    case 1:{ //fzzz
	      l1 = l2 = 0; l3 = 3;
	      break;
	    }
	    case 2:{ //fxzz
	      l1 = 1; l2 = 0; l3 = 2;
	      break;
	    }
	    case 3:{ //fyzz
	      l1 = 0; l2 = 1; l3 = 2;
	      break;
	    }
	    case 4:{ //fxxz
	      l1 = 2; l2 = 0; l3 = 1;
	      break;
	    }
	    case 5:{ //fyyz
	      l1 = 0; l2 = 2; l3 = 1;
	      break;
	    }
	    case 6:{ //fxxx
	      l1 = 3; l2 = l3 = 0;
	      break;
	    }
	    case 7:{ //fyyy
	      l1 = l3 = 0; l2 = 3;
	      break;
	    }
	    case 8:{ // fxyz
	      l1 = l2 = l3 = 1;
	      break;
	    }
	    case 9:{ // fxxy
	      l1 = 2; l2 = 1; l3 = 0;
	      break;
	    }
	    case 10:{ // fyyx
	      l1 = 1; l2 = 2; l3 = 0;
	      break;
	    }
	    }
	    break;
	  }
	  default:{
	    l1 = l2 = l3 = 0;
	  }
	  }
	  position = std::string::npos; // Exit loop
	}
      } else { // Move on to the next line
      	std::getline(input, line);
    }
    }
  } else {
    throw(Error("READBF", "Unable to open basis file."));
  }
  closeFile();
  BF b(c, l1, l2, l3, e);
  return b;
}

// Return a vector of the number of basis functions in each shell
// where the length of the vector is the number of shells
Vector BasisReader::readShells(int q)
{
   Vector shells(5); // Can't cope with higher than f functions, so 5 is max!	  
  // Open the file
  openFile(q);
  
  std::string delim = " ,";

  // Check it's open
  if (input.is_open()){
    // Parse
    std::string aname = getAtomName(q);
    aname += delim;
    std::string line, temp;
    int lmult;
    std::getline(input, line);
	int counter = 0;
    while(!input.eof()){
      if (line.find(aname) != std::string::npos){
	// Count the number of c lines that follow
	temp = line.at(0);
	int nbfs = 0;

	//Work out what shell type it is -
	// note that this calculates no. of cartesian gaussians
	if (temp == "s") { lmult = 1; }
	else if (temp == "sp") { lmult = 4; }
	else if (temp == "p") { lmult = 3; }
	else if (temp == "d") { lmult = 6; }
	else if (temp == "f") { lmult = 10; }

	std::getline(input, line);
	while (line.at(0) == 'c'){
	  nbfs += lmult;
	  std::getline(input, line);
	}
	shells[counter] = nbfs;
	counter++;
      } else { // Get next line
	std::getline(input, line);
      }
    }
  } else {
    throw(Error("IOERR", "Could not open basis file."));
  }

  // resize the return vector
  shells.resizeCopy(counter);
  
  // Close the file
  closeFile();
  
  return shells;
}

// Return a vector with the l-angular momentum quantum numbers of each
// shell belonging to atom q.
Vector BasisReader::readLnums(int q)
{
	Vector lnums(5); // See readShells
	  // Open the file
  openFile(q);
  
  std::string delim = " ,";

  // Check it's open
  if (input.is_open()){
    // Parse
    std::string aname = getAtomName(q);
    aname += delim;
    std::string line, temp;
    int lmult;
    std::getline(input, line);
	int counter = 0;
    while(!input.eof()){
      if (line.find(aname) != std::string::npos){
	// Count the number of c lines that follow
	temp = line.at(0);

	//Work out what shell type it is 
	if (temp == "s") { lmult = 0; }
	else if (temp == "sp") { lmult = 1; }
	else if (temp == "p") { lmult = 1; }
	else if (temp == "d") { lmult = 2; }
	else if (temp == "f") { lmult = 3; }

	lnums[counter] = lmult;
	counter++;
	std::getline(input, line);
      } else { // Get next line
	std::getline(input, line);
      }
    }
  } else {
    throw(Error("IOERR", "Could not open basis file."));
  }

  // resize the return vector
  lnums.resizeCopy(counter);
  
  // Close the file
  closeFile();
  
  return lnums;
}
	
	
// General purpose functions

// Return the atomic mass of an atom with atomic number q
double getAtomMass(int q)
{
	// Array of masses in atomic units for Hydrogen to Meitnerium
	// Ds and all the unun- types don't have reliable masses!
 	double masses[109] = { 1.0079, 4.0026, 6.941, 9.0122, 10.0811, 
 		12.0107, 14.0067, 15.9994, 18.9984, 20.1797, 22.9897,
 		24.305, 26.9815, 28.0855, 30.9738, 32.065, 35.453, 39.948,
 		39.0983, 40.078, 44.9559, 47.867, 50.9415, 51.9961, 54.938,
 		55.845, 58.9332, 58.6934, 63.546, 65.39, 69.723, 72.64,
 		74.9216, 78.96, 79.904, 83.8, 85.4678, 87.62, 88.9059, 91.224,
 		92.9064, 95.94, 98.0, 101.07, 102.9055, 106.42, 107.8682,
 		112.411, 114.818, 118.71, 121.76, 127.6, 126.9045, 131.293,
 		132.9055, 137.327, 138.9055, 140.116, 140.9077, 144.24, 145.0,
 		150.36, 151.964, 157.25, 158.9253, 162.5, 164.9303, 167.259,
 		168.9342, 173.04, 174.967, 178.49, 180.9479, 183.84, 186.207,
 		190.23, 192.217, 195.078, 196.9665, 200.59, 204.3833, 207.2,
 		208.9804, 209.0, 210.0, 222.0, 223.0, 226.0, 227.0, 232.0381,
 		231.0359, 138.0289, 237.0, 244.0, 243.0, 247.0, 247.0, 251.0,
 		252.0, 257.0, 258.0, 259.0, 262.0, 261.0, 262.0, 266.0, 264.0,
 		277.0, 268.0 }; 
 	return masses[q-1];
}

// Return the text version of atom with atomic number q -
// all caps is used for ease of parsing
// e.g., 6 -> C,   20 -> CA, etc.
std::string getAtomName(int q)
{
	std::string names[109] = {"H", "HE", "LI", "BE", "B", "C", "N",
		"O", "F", "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "AR",
		"K", "CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU",
		"ZN", "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR",
		"NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB",
		"TE", "I", "XE", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM",
		"EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA",
		"W", "RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO",
		"AT", "RN", "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM",
		"CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR", "RF", "DB", "SG",
		"BH", "HS", "MT" };
	return names[q-1];
}

// Return the atomic number of an atom with text n
// the reverse of getAtomName
int getAtomCharge(const std::string& n)
{
	// Make it upper case
	std::string name = n;
	std::transform(name.begin(), name.end(), name.begin(), ::toupper);
	int q = 0;
	bool found = false;
	while (!found && q < 109 ){
		  if(getAtomName(q) == name){
		  	found = true;
		  }
		  q++;
    }
    return q+1;
}

// Get the text name for a shell of angular momentum l
// e.g. l = 0 -> s,  l = 3 -> f
std::string getShellName(int l)
{
	std::string shell;
	switch(l) {
		case 0: {
			shell = "s";
			break;
		}
		case 1: {
			shell = "p";
			break;
		}	
		case 2: {
			shell = "d";
			break;
		}
		case 3: {
			shell = "f";
			break;
		}
		case 4: {
			shell = "g";
			break;
		}
		case 5: {
			shell = "h";
			break;
		}
		default: {
			shell = "N";
			break;
		}
	}		
	return shell;
}	