#include <chem/mcmm.h>
#include <chem/molecule.h>
#include <chem/mopac.h>
#include <chem/utils.h>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>

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
