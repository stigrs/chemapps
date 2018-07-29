////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2018 Stig Rune Sellevag. All rights reserved.
//
// This code is licensed under the MIT License (MIT).
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
////////////////////////////////////////////////////////////////////////////////

#include <chem/collision.h>
#include <chem/molecule.h>
#include <chem/thermochem.h>
#include <chem/thermodata.h>
#include <chem/troe.h>
#include <chem/whitten_rabino.h>
#include <srs/datum.h>
#include <srs/utils.h>
#include <boost/program_options.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

//
// Program for Troe factorization of low-pressure limiting
// rate coefficients.
//
int main(int argc, char* argv[])
{
    namespace po = boost::program_options;

    po::options_description options("Allowed options");
    // clang-format off
    options.add_options()
        ("help,h", "display help message")
        ("file,f", po::value<std::string>(), "input file");
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    std::string input_file;

    if (vm.find("help") != vm.end()) {
        std::cout << options << '\n';
        return 0;
    }
    if (vm.find("file") != vm.end()) {
        input_file = vm["file"].as<std::string>();
    }
    else {
        std::cerr << options << '\n';
        return 1;
    }

    try {
        std::ifstream from;
        srs::fopen(from, input_file);

        Molecule mol(from);
        Thermodata tpdata(from);
        Collision coll(from);
        Troe troe(from, mol);

        srs::Format<char> line;
        line.width(28).fill('-');

        std::cout << "Troe Factorization Analysis:\n"
                  << line('-') << '\n'
                  << "Abbreviations:\n"
                  << " k0^SC  - strong-collision low-pressure limiting rate"
                  << " coefficient\n"
                  << " Z_LJ   - Lennard-Jones collision frequency\n"
                  << " Q_vib  - vibrational partition function\n"
                  << " F_anh  - anharmonicity factor\n"
                  << " F_e    - energy dependence factor\n"
                  << " F_rot  - rotational factor\n"
                  << " F_free - free internal rotation factor\n"
                  << " F_hind - hindered internal rotation factor\n\n"
                  << "k0^SC and Z_LJ are given in cm^3 molecule^-1 s^-1\n\n";

        double e0  = troe.get_energy_barrier();
        double zpe = mol.get_vib().zero_point_energy();
        double rho = wr::vibr_density_states(mol, e0);
        double wra = wr::a_corr(mol, e0);

        std::cout << "Zero-point vibrational energy:    " << zpe << " cm^-1\n"
                  << "Vibrational density of states:    " << rho
                  << " (kJ/mol)^-1\n"
                  << "Whitten-Rabinovitch A correction: " << wra << "\n\n";

        auto temp = tpdata.get_temperature();
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}
