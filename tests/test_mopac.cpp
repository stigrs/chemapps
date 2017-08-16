#include <chem/datum.h>
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
        chem::fopen(from, "test_mopac.inp");
        Molecule mol(from);
        Mopac mop(from);
        mop.run(mol);

        const double heat_ans = 112.00281 * datum::cal_to_J;
        double heat           = mop.get_heat_of_formation();

        chem::Assert(std::abs(heat - heat_ans) < 1.0e-12,
                     std::runtime_error("bad heat of formation"));
    }
    catch (std::exception& e) {
        std::cerr << argv[0] << ": " << e.what() << '\n';
        return 1;
    }
}
