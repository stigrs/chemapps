#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <exception>

#include <chem/element.h>
#include <chem/utils.h>
#include <chem/molecule.h>
#include <armadillo>


int main(int argc, char* argv[])
{
    try {
        arma::vec3 rotc_ans = {127.63201, 24.89071, 24.02767};

        std::ifstream from;
        chem::fopen(from, "test_molrot.inp");

        Molecule mol(from, std::cout, "Molecule");
        arma::vec3 rotc = mol.get_rot()->constants();

        for (int i = 0; i < rotc.size(); ++i) {
            chem::Assert(std::abs(rotc(i) - rotc_ans(i)) < 1.0e-4, 
                         std::runtime_error("bad rot const"));
        }
    }
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
        return 1;
    }
}
