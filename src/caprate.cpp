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

#include <chem/traits.h>
#include <chem/molecule.h>
#include <chem/thermochem.h>
#include <numlib/matrix.h>
#include <numlib/constants.h>
#include <numlib/math.h>
#include <stdutils/stdutils.h>
#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <memory>

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

Numlib::Grid e_grid;                   // energy grid
Numlib::Grid j_grid;                   // angular momentum grid
Numlib::Grid t_grid;                   // temperature grid
Numlib::Mat<double> nej;               // N(E,J) data
std::unique_ptr<Chem::Molecule> frag1; // input data on fragment 1
std::unique_ptr<Chem::Molecule> frag2; // input data on fragment 2

//-----------------------------------------------------------------------------

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
    Stdutils::fopen(from, input_file);

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

    frag1 = std::make_unique<Chem::Molecule>(from, "Fragment1", true);
    frag1->rot_analysis();
    std::cout << '\n';
    frag2 = std::make_unique<Chem::Molecule>(from, "Fragment2", true);
    frag2->rot_analysis();
}

void read_nej(const std::string& input_file)
{
    std::ifstream from;
    Stdutils::fopen(from, input_file);

    nej.resize(e_grid.size(), j_grid.size());

    double ee = 0.0;
    double jj = 0.0;
    double n = 0.0;
    for (Index j = 0; j < j_grid.size(); ++j) {
        for (Index e = 0; e < e_grid.size(); ++e) {
            from >> ee >> jj >> n;
            if (!from) {
                throw Caprate_error("cannot read " + input_file);
            }
            if (ee != e_grid[e]) {
                throw Caprate_error(
                    std::to_string(e) + "-th E has bad value: " +
                    std::to_string(ee) + ", " + std::to_string(e_grid[e]));
            }
            nej(e, j) = n;
        }
        if (jj != j_grid[j]) {
            throw Caprate_error(std::to_string(j) +
                                "-th J has bad value: " + std::to_string(jj) +
                                ", " + std::to_string(j_grid[j]));
        }
    }
}

// Integrate N(E,J) to obtain capture rate coefficient.
//
// The fragment partition function is given by the relative translational
// partition function of the fragments and the rotational partition functions
// of each fragment.
//
void integrate()
{
    namespace Pc = Numlib::Constants;

    std::cout << "\nAbbreviations:\n"
              << " Q_frag(T) - partition function of fragments (cm**-3).\n"
              << " k_cap(T)  - capture rate coefficient "
              << "(cm**3 molecule**-1 s**-1).\n\n"
              << "T/K\tQ_frag(T)\tk_cap(T)\n"
              << std::setw(35) << std::setfill('-') << '-' << std::setfill(' ')
              << '\n';

    double estart = e_grid.start() / Pc::au_to_icm;
    double emax = e_grid.size() / Pc::au_to_icm;
    double estep = e_grid.step() / Pc::au_to_icm;
    double jstep = j_grid.step();

    e_grid.set(estart, emax, estep); // atomic units

    double tt;
    double twoj;
    double qfrag;
    double kcap;

    double red_mass = frag1->tot_mass() * frag2->tot_mass() /
                      (frag1->tot_mass() + frag2->tot_mass());

    for (Index t = 0; t < t_grid.size(); ++t) {
        kcap = 0.0;
        tt = t_grid[t] / Pc::au_to_K;
        for (Index e = 0; e < e_grid.size(); ++e) {
            for (Index j = 0; j < j_grid.size(); ++j) {
                twoj = 2.0 * j_grid[j] + 1.0;
                kcap += twoj * nej(e, j) * std::exp(-e_grid[e] / tt);
            }
        }
        qfrag = Chem::qtrans(red_mass, t_grid[t]) * 1.0e-6; // cm**-3
        qfrag *= Chem::qrot(*frag1, t_grid[t], true);
        qfrag *= Chem::qrot(*frag2, t_grid[t], true);
        kcap *= estep * jstep / (2.0 * Pc::pi * Pc::h_bar * qfrag / Pc::E_h);

        std::cout << t_grid[t] << '\t' << qfrag << '\t' << kcap << '\n';
    }
}

