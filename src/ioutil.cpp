/*
 *
 *   PURPOSE: Implements ioutil.hpp, a collection of i/o utility functions.
 *
 *   DATE           AUTHOR           CHANGES
 *   ==================================================================
 *   30/08/15       Robert Shaw      Original code.
 *
 */

// Includes
#include "ioutil.hpp"
#include <algorithm>	
	
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
	int q = 1;
	bool found = false;
	while (!found && q < 109 ){
		  if(getAtomName(q) == name){
		  	found = true;
		  }
		  q++;
    }
    return q-1;
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
