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

#ifndef CHEM_GAMSS_H
#define CHEM_GAMSS_H

#include <chem/molecule.h>
#include <chem/conformer.h>
#include <iostream>
#include <string>
#include <vector>
#include <random>

namespace Chem {

// Class providing molecular structure search using a genetic algorithm.
//
template <class Pot>
class Gamss {
public:
    Gamss(std::istream& from, std::ostream& to = std::cout);

    // Run solver.
    void solve(std::ostream& to = std::cout);

private:
    // Initialize population.
    void init_population(std::ostream& to = std::cout);

    // Generate a new random conformer.
    void gen_rand_conformer(Molecule& m);

    // Select random dihedral angle.
    std::vector<int> select_rand_dihedral(const Molecule& m);

    // Check if geometry of random conformer is sensible.
    bool geom_sensible(const Molecule& m) const;

    // Check if geometry is blacklisted.
    bool is_blacklisted(const Numlib::Mat<double>& xyz) const;

    Molecule mol; // molecule
    Pot pot;      // potential function

    double dist_min; // smallest atom-atom distance permitted
    double dist_max; // largest bond distance permitted

    double energy_min; // smallest energy permitted
    double energy_max; // largest energy permitted

    double rmsd_tol_uniq; // geometry tolerance

    int pop_size;     // population size
    int max_iter;     // max number of iterations
    int max_mut_tors; // max number of torsional mutations
    int seed;         // random number generator seed

    std::vector<Conformer> population; // population of optimized structures
    std::vector<Conformer> blacklist;  // population of blacklisted structures

    std::mt19937_64 mt; // random number engine
};

} // namespace Chem

#endif // CHEM_GAMSS_H

