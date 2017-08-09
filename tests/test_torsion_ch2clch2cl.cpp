#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include <exception>

#include <chem/utils.h>
#include <chem/molecule.h>


int main(int /* argc */, char* argv[])
{
    try {
        std::ifstream from;
        chem::fopen(from, "test_torsion_ch2clch2cl.inp");

        Molecule mol(from, std::cout, "Molecule");
        double rmi = mol.get_tor()->red_moment_of_inertia();

        const double rmi_ans = 58.76991427; // Chuang and Truhlar (2000)
        chem::Assert(std::abs(rmi - rmi_ans) < 9.0e-2,
                     std::runtime_error("bad rmi"));
    }
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
        return 1;
    }
}
