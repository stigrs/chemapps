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
#include <iomanip>
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
			      << "Energy barrier towards reaction:  " << e0 << " cm^-1\n"
                  << "Vibrational density of states:    " << rho
                  << " (kJ/mol)^-1\n"
                  << "Whitten-Rabinovitch A correction: " << wra << "\n\n";

        line.width(75).fill('-');

        // clang-format off
        std::cout << std::setiosflags(std::ios_base::left) << line('-') << '\n'
                  << std::setw(8)  << "T/K" 
			      << std::setw(11) << "k0^SC/[M]"
                  << std::setw(10) << "Z_LJ/[M]" 
			      << std::setw(8)  << "Q_vib"
                  << std::setw(8)  << "F_anh" 
			      << std::setw(8)  << "F_e"
                  << std::setw(8)  << "F_rot" 
			      << std::setw(8)  << "F_free"
                  << std::setw(8)  << "F_hind" 
			      << '\n' << line('-') << '\n'
                  << std::resetiosflags(std::ios_base::left);
        // clang-format on

        srs::Format<double> gen65;
        srs::Format<double> sci82;
        srs::Format<double> sci92;

        gen65.width(6).precision(5);
        sci82.scientific().width(8).precision(2);
        sci92.scientific().width(9).precision(2);

        auto temp = tpdata.get_temperature();

        for (auto t : temp) {
            double kT     = datum::R * 1.0e-3 * t;
            double z_lj   = coll.lj_coll_freq(t);
            double q_vib  = chem::qvib(mol, t, "V=0");
            double f_anh  = troe.f_anharm();
            double f_e    = troe.f_energy(t);
            double f_rot  = troe.f_rotation(t);
            double f_free = troe.f_free_rotor(t);
            double f_hind = troe.f_hind_rotor(t);

            double k0 = z_lj * (rho * kT / q_vib) * f_anh * f_e * f_rot * f_free
                        * f_hind * std::exp(-e0 * datum::icm_to_kJ / kT);

            // clang-format off
            std::cout << gen65(t)      << "  " 
				      << sci92(k0)     << "  " 
				      << sci82(z_lj)   << "  " 
				      << gen65(q_vib)  << "  " 
				      << gen65(f_anh)  << "  "
                      << gen65(f_e)    << "  " 
				      << gen65(f_rot)  << "  "
                      << gen65(f_free) << "  " 
				      << gen65(f_hind) << '\n';
            // clang-format on
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}
