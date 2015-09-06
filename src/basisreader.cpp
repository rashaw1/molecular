/*
 *
 *   PURPOSE: Implements basisreader.hpp, a class for reading basis files.
 *
 *   DATE           AUTHOR           CHANGES
 *   ==================================================================
 *   30/08/15       Robert Shaw      Original code.
 *   04/09/15       Robert Shaw      Now indexes primitives.
 */
 
 #include "basisreader.hpp"
 #include "ioutil.hpp"
 #include "vector.hpp"
 #include "bf.hpp"
 #include "error.hpp"
 #include <algorithm>
#include <iostream>
 
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
  Vector c; Vector e; Vector ids; // Coeffs, exps, and prim ids

  // Check it's open
  if (input.is_open()){
    // Parse
    std::string aname = getAtomName(q);
    aname += " "; aname += delim;
    std::size_t position;
    int bfcount = -1;
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
	long nexps = std::count(line.begin(), line.end(), ',') - 1;
	// Copy in
	Vector tempexps((int) nexps);
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
	  std::size_t p = temp2.find('.');
	  int start = std::stoi(temp2.substr(0, p));
	  int end = std::stoi(temp2.substr(p+1, temp2.length()));
	  
	  // Now resize e and ids, and copy in
	  e.resize(end - start + 1);
	  ids.resize(end - start + 1);
	  for (int j = 0; j < end - start + 1; j++){
	    e[j] = tempexps(start+j-1); // Take account of zero-indexing
	    ids[j] = start*lmult + j - lmult;
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
	    case 1: { //pz 
	      l3 = 1; l2 = l1 = 0;
	      break;
	    }
	    case 2: { //py
	      l1 = l3 = 0; l2 = 1;
	      for (int index = 0; index < ids.size(); index++) { ids[index] += e.size(); }
	      break;
	    }
	    case 3: { //px
	      l3 = l2 = 0; l1 = 1;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 2*e.size(); }
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
	    case 2:{ // pz
	      l3 = 1; l2 = l1 = 0;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += e.size(); }
	      break;
	    }
	    case 3:{ // py
	      l1 = l3 = 0; l2 = 1;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 2*e.size(); }
	      break;
	    }
	    case 4:{ //px
	      l3 = l2 = 0; l1 = 1;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 3*e.size(); }
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
	    case 2:{ // dyz
	      l3 = l2 = 1; l1 = 0;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += e.size(); }
	      break;
	    }
	    case 3:{ // dxz
	      l1 = l3 = 1; l2 = 0;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 2*e.size(); }
	      break;
	    }
	    case 4:{ // dyy
	      l2 = 2; l3 = l1 = 0;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 3*e.size(); }
	      break;
	    }
	    case 5:{ // dxy 
	      l3 = 0; l2 = l1 = 1;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 4*e.size(); }
	      break;
	    }
	    case 6:{ // dxx
	      l1 = 2; l2 = l3 = 0;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 5*e.size(); }
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
	    case 2:{ //fyzz
	      l1 = 0; l2 = 1; l3 = 2;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += e.size(); }
	      break;
	    }
	    case 3:{ //fxzz
	      l1 = 1; l2 = 0; l3 = 2;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 2*e.size(); }
	      break;
	    }
	    case 4:{ //fyyz
	      l1 = 0; l2 = 2; l3 = 1;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 3*e.size(); }
	      break;
	    }
	    case 5:{ //fxyz
	      l1 = 1; l2 = 1; l3 = 1;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 4*e.size(); }
	      break;
	    }
	    case 6:{ //fxxz
	      l1 = 2; l2 = 0; l3 = 1;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 5*e.size(); }
	      break;
	    }
	    case 7:{ //fyyy
	      l1 = l3 = 0; l2 = 3;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 6*e.size(); }
	      break;
	    }
	    case 8:{ // fxyy
	      l1 = 1; l2 = 2; l3 = 0;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 7*e.size(); }
	      break;
	    }
	    case 9:{ // fxxy
	      l1 = 2; l2 = 1; l3 = 0;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 8*e.size(); }
	      break;
	    }
	    case 10:{ // fxxx
	      l1 = 3; l2 = 0; l3 = 0;
	      for (int index = 0; index< ids.size(); index++) { ids[index] += 9*e.size(); }
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
  BF b(c, l1, l2, l3, e, ids);
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
  int counter = 0;
  // Check it's open
  if (input.is_open()){
    // Parse
    std::string aname = getAtomName(q);
    aname += delim;
    std::string line, temp;
    int lmult;
    std::getline(input, line);

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
  int counter = 0;
  // Check it's open
  if (input.is_open()){
    // Parse
    std::string aname = getAtomName(q);
    aname += delim;
    std::string line, temp;
    int lmult;
    std::getline(input, line);

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

// Now do these for a complete set of qs
Vector BasisReader::readShells(Vector& qs)
{
  Vector shells;
  Vector temp;
  shells = readShells(qs(0));
  for (int i = 1; i < qs.size(); i++){
    int s = shells.size();
    temp = readShells(qs(i));
    int t = temp.size();
    shells.resizeCopy(s+t);
    for (int j = 0; j < t; j++) {
      shells[j+s] = temp[j];
    }
  }
  return shells;
}

Vector BasisReader::readLnums(Vector& qs)
{
  Vector lnums;
  Vector temp;
  lnums = readLnums(qs(0));
  for (int i = 1; i < qs.size(); i++){
    int l = lnums.size();
    temp = readLnums(qs(i));
    int t = temp.size();
    lnums.resizeCopy(l+t);
    for (int j = 0; j < t; j++){
      lnums[j+l] = temp[j];
    }
  }
  return lnums;
}
