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

#ifndef CHEM_MCMM_H
#define CHEM_MCMM_H

#include <chem/conformer.h>
#include <chem/molecule.h>
#include <numlib/matrix.h>
#include <numlib/math.h>
#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace Chem {

// Class providing Monte Carlo Multiple Minima (MCMM) solver.
//
template <class Pot>
class Mcmm {
public:
    Mcmm(std::istream& from,
         Molecule& mol_,
         const std::string& key = "Mcmm",
         bool verbose_ = false);

    // MCMM solver.
    void solve(std::ostream& to = std::cout);

    // Get global minimum values.
    double get_global_min_energy();
    Numlib::Mat<double> get_global_min_xyz();

private:
    // Check if MCMM solver is finished.
    bool check_exit() const;

    // Check acceptance of energy.
    bool accept_energy(double enew);

    // Check if random conformer is ok.
    bool accept_geom_dist(const Molecule& m) const;

    // Check if current conformer is a duplicate.
    bool duplicate(const Molecule& m) const;

    // Generate a new molecular conformer.
    void new_conformer();

    // Update MCMM solver.
    void update();

    // Function for selecting starting geometry by using the uniform
    // usage scheme.
    void uniform_usage(Numlib::Mat<double>& xnew);

    // Generate a new random conformer.
    void gen_rand_conformer(Molecule& m);

    // Save conformer (local energy minimum).
    void save_conformer(const Molecule& m);

    // Sort conformers in ascending order and cut to nminima.
    void sort_conformers();

    // Select random dihedral angle.
    std::vector<int> select_rand_dihedral(const Molecule& m);

    Molecule& mol; // molecule
    Pot pot;       // potential function

    double xtol; // absolute error in geometry
    double etol; // absolute error in energy
    double emin; // lowest energy permitted
    double emax; // highest energy permitted
    double rmin; // smallest bond distance permitted
    double temp; // temperature

    unsigned maxiter;   // maximum number of iterations (k)
    unsigned miniter;   // minimum number of iterations (k)
    unsigned maxreject; // maximum number of consecutive rejected trials
    unsigned nminima;   // number of local minima stored

    int seed; // random number generator seed

    unsigned kiter;   // iteration parameter
    unsigned nreject; // iterator for number of consecutive rejected trials
    unsigned naccept; // iterator for number of accepted trials

    Numlib::Mat<double> xcurr;   // current geometry
    Numlib::Mat<double> xglobal; // geometry of global minimum

    double ecurr;                // current energy
    std::vector<double> eglobal; // energy of global minimum

    std::vector<Conformer> conformers; // array with local energy minima

    bool verbose;
    bool global_min_found = false;

    std::mt19937_64 mt; // random number engine
};

template <class Pot>
inline double Mcmm<Pot>::get_global_min_energy()
{
    if (!global_min_found) {
        solve();
    }
    return *std::min_element(eglobal.begin(), eglobal.end());
}

template <class Pot>
inline Numlib::Mat<double> Mcmm<Pot>::get_global_min_xyz()
{
    if (!global_min_found) {
        solve();
    }
    return xglobal;
}
template <class Pot>
inline void Mcmm<Pot>::gen_rand_conformer(Molecule& m)
{
    // Select a random dihedral angle:
    std::vector<int> moiety = select_rand_dihedral(m);

    // Apply random variation to dihedral angle:
    std::uniform_real_distribution<> rnd_uni_real(-180.0, 180.0);
    double delta = rnd_uni_real(mt);
    m.int_coord().rotate_moiety(moiety, delta);
}

template <class Pot>
void Mcmm<Pot>::save_conformer(const Molecule& m)
{
    Conformer c(m.elec_energy(), m.cart_coord());
    conformers.push_back(c);
}

} // namespace Chem

#endif // CHEM_MCMM_H
