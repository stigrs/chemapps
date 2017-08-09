#include <iostream>
#include <fstream>
#include <stdexcept>
#include <exception>
#include <chem/mopac.h>
#include <chem/utils.h>


int main(int /*argc */, char* argv[])
{
    try {
        std::ifstream from;
        chem::fopen(from, "test_mopac.inp");
        Mopac mop(from, "Mopac");
    }
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
        return 1;
    }
}
