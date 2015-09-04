#include <iostream>
#include <fstream>
#include <iomanip>
#include "logger.hpp"
#include "molecule.hpp"
#include "error.hpp"
#include "integrals.hpp"

int main (int argc, char* argv[])
{
  std::ifstream input("inputfile.mol");
  if (input.is_open()) { std::cout << "input opened\n"; }
  std::ofstream output("out.out");
  if (output.is_open()) { std::cout << "output opened\n"; }
  Logger log(input, output, std::cout);
  std::cout << "log made\n";
  try {
  Molecule mol(log, 1);
  std::cout << "molecule made\n";

  log.print(mol, true);
  log.print(log.getBasis(), true);
  //log.print(mol.getBF(7, 6));
  //log.print(mol.getBF(7, 6).getPBF(0));
  log.result("This is a result.");
  std::cout << "All fine\n";
  IntegralEngine integral(mol);
  std::cout << integral.getOverlap(41, 35) << "  " << integral.getKinetic(41, 35) << "\n";
  log.print("\n");
  log.print(mol.getAtom(2).getBF(35));
  log.print("\n");
  log.print(mol.getAtom(2).getBF(29));
  log.print("\n");
  for(int i = 0; i < 4; i++){
    log.print("\n");
    log.print(mol.getAtom(2).getBF(29).getPBF(i));
  }
  Vector ests;
  ests = (1.0/(1024.0*1024.0))*integral.getEstimates();
  log.print(ests);
  log.print("\n");
  log.finalise();
  } catch (Error e) {
    log.error(e);
  }
  input.close();
  std::cout << "input closed\n";
  output.close();
  std::cout << "output closed\n";
  return 0;
}
