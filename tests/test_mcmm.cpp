#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <exception>
#include <chem/molecule.h>
#include <chem/mopac.h>
#include <chem/mcmm.h>
#include <chem/utils.h>


int main(int /*argc */, char* argv[])
{
    try {
        std::ifstream from;
        chem::fopen(from, "test_mcmm.inp");

        Molecule mol;
        Mcmm<Mopac> mc(from, mol);
    }
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
        return 1;
    }
}
