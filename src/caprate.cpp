////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2017 Stig Rune Sellevag. All rights reserved.
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

#include <chem/mol_type.h>
#include <chem/molecule.h>
#include <chem/thermochem.h>
#include <srs/array.h>
#include <srs/datum.h>
#include <srs/math.h>
#include <srs/utils.h>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

//-----------------------------------------------------------------------------

// Error reporting:

struct Caprate_error : std::runtime_error {
    Caprate_error(std::string s) : std::runtime_error(s) {}
};

//-----------------------------------------------------------------------------

// Forward declarations:

void integrate();
void read_input(const std::string& input_file);
void read_nej(const std::string& input_file);

//-----------------------------------------------------------------------------

// Global declarations:

Grid e_grid;                // energy grid
Grid j_grid;                // angular momentum grid
Grid t_grid;                // temperature grid
Molecule frag1;             // input data on fragment 1
Molecule frag2;             // input data on fragment 2
srs::Array<double, 2> nej;  // N(E,J) data

//-----------------------------------------------------------------------------

//
// Program for obtaining capture rate coefficient of the association reaction
//
//         A + B ---> AB
//
// by integration of N(E,J) obtained from a VRC-TST calculation. The
// translational, rotational, and vibrational degrees of freedom are assumed
// to be divided into sets of transitional and conserved degrees of freedom.
//
int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cerr << "usage: " << argv[0] << " input_file nej_file\n";
        return 1;
    }
    try {
        read_input(argv[1]);
        read_nej(argv[2]);
        integrate();
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

//-----------------------------------------------------------------------------

void read_input(const std::string& input_file)
{
    std::ifstream from;
    srs::fopen(from, input_file);

    e_grid.set(from, "EnergyGrid");
    j_grid.set(from, "AngMomGrid");
    t_grid.set(from, "TemprsGrid");

    std::cout << '\n'
              << "Specification of E grid (cm**-1):\n"
              << e_grid << '\n'
              << "Specification of J grid (au):\n"
              << j_grid << '\n'
              << "Specification of T grid (K):\n"
              << t_grid << '\n';

    frag1.init(from, std::cout, "Fragment1", true);
    frag1.get_rot().analysis();
    std::cout << '\n';
    frag2.init(from, std::cout, "Fragment2", true);
    frag2.get_rot().analysis();
}

void read_nej(const std::string& input_file)
{
    std::ifstream from;
    srs::fopen(from, input_file);

    nej.resize(e_grid.size(), j_grid.size());

    double ee = 0.0;
    double jj = 0.0;
    double n  = 0.0;
    for (srs::size_t j = 0; j < j_grid.size(); ++j) {
        for (srs::size_t e = 0; e < e_grid.size(); ++e) {
            from >> ee >> jj >> n;
            if (!from) {
                throw Caprate_error("cannot read " + input_file);
            }
            if (ee != e_grid[e]) {
                throw Caprate_error(srs::to_string(e) + "-th E has bad value: "
                                    + srs::to_string(ee) + ", "
                                    + srs::to_string(e_grid[e]));
            }
            nej(e, j) = n;
        }
        if (jj != j_grid[j]) {
            throw Caprate_error(srs::to_string(j)
                                + "-th J has bad value: " + srs::to_string(jj)
                                + ", " + srs::to_string(j_grid[j]));
        }
    }
}

//
// Integrate N(E,J) to obtain capture rate coefficient.
//
// The fragment partition function is given by the relative translational
// partition function of the fragments and the rotational partition functions
// of each fragment.
//
void integrate()
{
    std::cout << "\nAbbreviations:\n"
              << " Q_frag(T) - partition function of fragments (cm**-3).\n"
              << " k_cap(T)  - capture rate coefficient "
              << "(cm**3 molecule**-1 s**-1).\n\n"
              << "T/K\tQ_frag(T)\tk_cap(T)\n"
              << std::setw(35) << std::setfill('-') << '-' << std::setfill(' ')
              << '\n';

    double estart = e_grid.start() / datum::au_to_icm;
    double emax   = e_grid.size() / datum::au_to_icm;
    double estep  = e_grid.step() / datum::au_to_icm;
    double jstep  = j_grid.step();

    e_grid.set(estart, emax, estep);  // atomic units

    double tt;
    double twoj;
    double qfrag;
    double kcap;

    double red_mass = frag1.tot_mass() * frag2.tot_mass()
                      / (frag1.tot_mass() + frag2.tot_mass());

    for (srs::size_t t = 0; t < t_grid.size(); ++t) {
        kcap = 0.0;
        tt   = t_grid[t] / datum::au_to_K;
        for (srs::size_t e = 0; e < e_grid.size(); ++e) {
            for (srs::size_t j = 0; j < j_grid.size(); ++j) {
                twoj = 2.0 * j_grid[j] + 1.0;
                kcap += twoj * nej(e, j) * std::exp(-e_grid[e] / tt);
            }
        }
        qfrag = chem::qtrans(red_mass, t_grid[t]) * 1.0e-6;  // cm**-3
        qfrag *= chem::qrot(frag1, t_grid[t], true);
        qfrag *= chem::qrot(frag2, t_grid[t], true);
        kcap *= estep * jstep
                / (2.0 * datum::pi * datum::h_bar * qfrag / datum::E_h);

        std::cout << t_grid[t] << '\t' << qfrag << '\t' << kcap << '\n';
    }
}
