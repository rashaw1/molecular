/*
 *
 *   PURPOSE: Implements filereader.hpp, a class for reading input files.
 *
 *   DATE           AUTHOR           CHANGES
 *   ==================================================================
 *   30/08/15       Robert Shaw      Original code.
 *
 */
 
 #include "filereader.hpp"
 #include "ioutil.hpp"
 #include "error.hpp"
 #include <algorithm>
#include <iostream> 

 // Class FileReader implementation

// Destructor
FileReader::~FileReader()
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
  else if (t == "thrint") { rval = 7; }
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
  thrint = 1e-12;

  // Read line by line and parse
  std::string line, token;
  size_t pos;
  int linecount = 0;
  while (std::getline(input, line)) {
    // Erase any comments
    pos = line.find('!');
    if (pos != std::string::npos){
      line.erase(pos, line.length());
    }

    // Tokenise
    pos = line.find(',');
    if(pos != std::string::npos){
      token = line.substr(0, pos);

      // Match token to command
      switch(findToken(token)){
      case 1: { // Charge
	charge = std::stoi(line.substr(pos+1, line.length()));
	break;
      }
      case 2: { // Multiplicity
	multiplicity = std::stoi(line.substr(pos+1, line.length()));
	break;
      }
      case 3: { // MaxIter
	maxiter = std::stoi(line.substr(pos+1, line.length()));
	break;
      }
      case 4: { // Precision
	precision = std::stod(line.substr(pos+1, line.length()));
	break;
      }
      case 5: { // Basis
	line.erase(0, pos+1);
	// Get rid of excess spaces
	line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
	basis = line;
	break;
      }
      case 6: { // GeomStart
	geomstart = linecount + 1;
	geomend = geomstart-1;
	while(line != "geomend"){
	  std::getline(input, line);
	  geomend++;
	}
	break;
      }
      case 7: { // Thrint
	thrint = std::stod(line.substr(pos+1, line.length()));
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
    std::string line;
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
