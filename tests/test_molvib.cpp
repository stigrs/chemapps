#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <exception>

#include <chem/element.h>
#include <chem/utils.h>
#include <chem/datum.h>
#include <chem/molecule.h>
#include <armadillo>


int main(int /* argc */, char* argv[])
{
    try {
        double zpe_ans = 0.052023;

        std::ifstream from;
        chem::fopen(from, "test_molvib.inp");

        Molecule mol(from, std::cout, "Molecule");
        double zpe = mol.get_vib()->zero_point_energy() / datum::au2icm;

        chem::Assert(std::abs(zpe - zpe_ans) < 1.0e-6, 
                     std::runtime_error("bad zero-point energy"));
    }
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
        return 1;
    }
}
