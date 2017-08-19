///////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////

#ifndef CHEM_MCMM_H
#define CHEM_MCMM_H

#include <chem/conformer.h>
#include <chem/math.h>
#include <chem/molecule.h>
#include <armadillo>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// Error reporting:

struct Mcmm_error : std::runtime_error {
    Mcmm_error(std::string s) : std::runtime_error(s) {}
};

//
// Class providing Monte Carlo Multiple Minima (MCMM) solver.
//
template <class Pot>
class Mcmm {
public:
    Mcmm(std::istream& from,
         Molecule& mol_,
         const std::string& key = "Mcmm",
         bool verbose_          = false);

    ~Mcmm() {}

    // MCMM solver.
    void solve();

private:
    // Generate a new molecular conformer.
    void new_conformer();

    // Update MCMM solver.
    void update();

    // Check if MCMM solver is finished.
    bool check_exit() const;

    // Check if random conformer is ok.
    bool accept_geom_dist(const Molecule& m) const;

    // Function for selecting starting geometry by using the uniform
    // usage scheme.
    void uniform_usage(arma::mat& xnew);

    // Generate a new random conformer.
    void gen_rand_conformer(Molecule& m);

    // Select random dihedral angle.
    std::vector<int> select_rand_dihedral(const Molecule& m);

#if 0
    // Check acceptance of energy.
    bool accept_energy(double enew);

    // Sort local minima in descending order and cut to nminima.
    void sort_minima();

    // Save local minima.
    void save_conformer(double energy, const Molecule& mol);
#endif
    Molecule& mol;  // molecule
    Pot pot;        // potential function

    double xtol;  // absolute error in geometry
    double etol;  // absolute error in energy
    double emin;  // lowest energy permitted
    double emax;  // highest energy permitted
    double rmin;  // smallest bond distance permitted
    double temp;  // temperature

    int maxiter;    // maximum number of iterations (k)
    int maxreject;  // maximum number of consecutive rejected trials
    int nminima;    // number of local minima stored
    int seed;       // random number generator seed

    int kiter;    // iteration parameter
    int nreject;  // iterator for number of consecutive rejected trials
    int naccept;  // iterator for number of accepted trials

    arma::mat xcurr;    // current geometry
    arma::mat xglobal;  // geometry of global minimum

    double ecurr;    // current energy
    double eglobal;  // energy of global minimum

    std::vector<Conformer> conformers;  // array with local energy minima

    bool verbose;

    std::mt19937_64 mt;  // random number engine
};

template <class Pot>
inline void Mcmm<Pot>::gen_rand_conformer(Molecule& m)
{
    // Select a random dihedral angle:
    std::vector<int> moiety = select_rand_dihedral(m);

    // Apply random variation to dihedral angle:
    std::uniform_real_distribution<> rnd_uni_real(-180.0, 180.0);
    double delta = rnd_uni_real(mt);
    m.get_zmat()->rotate_moiety(moiety, delta);
}

template <class Pot>
inline bool Mcmm<Pot>::accept_geom_dist(const Molecule& m) const
{
    bool geom_ok = true;
    arma::mat dist_mat;
    chem::pdist_matrix(dist_mat, m.get_xyz());
    // std::cout << arma::nonzeros(dist_mat).min() << '\n';
    if (arma::nonzeros(dist_mat).min() < rmin) {  // avoid too close atoms
        geom_ok = false;
    }
    return geom_ok;
}

#endif  // CHEM_MCMM_H
