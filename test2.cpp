#include <iostream>
#include <fstream>
#include <iomanip>
#include "mathutil.hpp"
#include "vector.hpp"
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
  std::cout << mol.getAtom(2).getNShellPrims(1);
  std::cout << "\n";
  std::cout << mol.getAtom(0).getBF(0).getPBF(0).getNorm() << " ";
  std::cout << mol.getAtom(0).getBF(0).getPBF(1).getNorm() << " ";
  std::cout << mol.getAtom(0).getBF(0).getPBF(2).getNorm() << "\n";

  Vector c = mol.getAtom(0).getCoords();  
  Vector cd = mol.getAtom(1).getCoords();
  Vector d = mol.getAtom(2).getCoords();

  IntegralEngine integral(mol);
  std::cout << "Fine\n";
  std::cout << "Nuc Attract = " << integral.nucAttract(mol.getAtom(0).getBF(0).getPBF(2), mol.getAtom(2).getBF(0).getPBF(2), 
						       c, d, cd) << "\n";
  std::cout << integral.getOverlap(6, 6) << "  " << integral.getKinetic(6, 6) << "\n";
    log.print("\n");
  //log.print(mol.getAtom(2).getBF(28));
  //log.print("\n");
  //log.print(mol.getAtom(2).getBF(24));
  //log.print("\n");
  Vector ests;
  ests = (1.0/(1024.0*1024.0))*integral.getEstimates();
  log.print("\n");log.print(ests);log.print("\n");
  Matrix testing(3, 3, 0.5);
  for (int i = 0; i < 3; i++) { testing(i, i) = 1.0; }
  testing = integral.makeSpherical(1, 1, testing);
  testing.print();
  std::cout << "\n";
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
