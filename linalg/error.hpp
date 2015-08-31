/*
 *     PURPOSE: defines and implements a simple error class, containing
 *              an error code, and a message, and a way of printing
 *              to the primary output stream.
 * 
 *     DATE             AUTHOR                CHANGES
 *   =======================================================================
 *     14/08/15         Robert Shaw           Original code
 */

#ifndef ERRORHEADERDEF
#define ERRORHEADERDEF

// Includes
#include <string>
#include <map>

class Error
{
private:
  std::string code; // Shorthand code for error
  std::string msg; // Detailed error message
public:
  // Constructors
  Error(); // Default constructor - general, unspecific error
  Error(std::string c, std::string m); // Specify custom code and message
  // Look up the error message from pre-specified list errormap
  // using only the shorthand code:
  Error(std::string c, std::map<std::string, std::string> errormap);
  // Print to primary output stream, full = false gives just the code
  void print(bool full = true); 
  // Accessors
  std::string getCode() const { return code; }
  std::string getMsg() const { return msg; }
};

#endif
