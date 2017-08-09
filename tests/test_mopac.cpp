#include <iostream>
#include <fstream>
#include <stdexcept>
#include <exception>
#include <chem/mopac.h>
#include <chem/molecule.h>
#include <chem/utils.h>


int main(int /*argc */, char* argv[])
{
    try {
        std::ifstream from;
        chem::fopen(from, "test_mopac.inp");
        Molecule mol(from, std::cout, "Molecule");
        Mopac mop(from, "Mopac");
        mop.run(mol);
    }
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
        return 1;
    }
}
