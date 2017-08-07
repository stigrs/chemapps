#include <iostream>
#include <fstream>
#include <stdexcept>
#include <exception>

#include <chem/utils.h>
#include <chem/molecule.h>
#include <armadillo>


int main(int /* argc */, char* argv[])
{
    try {
        std::ifstream from;
        chem::fopen(from, "test_torsion_ch3oh.inp");

        Molecule mol(from, std::cout, "Molecule", true);
    }
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
        return 1;
    }
}
