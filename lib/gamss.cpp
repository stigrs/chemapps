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

#include <chem/gamss.h>
#include <chem/gaussian.h>
#include <chem/mopac.h>
#include <chem/impl/io_support.h>
#include <numlib/matrix.h>
#include <numlib/math.h>
#include <stdutils/stdutils.h>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <limits>

template <class Pot>
Chem::Gamss<Pot>::Gamss(std::istream& from, std::ostream& to) : mol(from, to)
{
    using namespace Stdutils;
    const std::string key = "Gamss";

    // Read input:

    dist_min = 0.5;
    dist_max = 2.5;

    energy_min = -std::numeric_limits<double>::max();
    energy_max = 0.0;

    rmsd_tol_uniq = 0.2;

    pop_size = 20;
    max_iter = 200;
    max_mut_tors = 2;
    seed = 0;

    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "dist_min", dist_min, dist_min);
        get_token_value(from, pos, "dist_max", dist_max, dist_max);
        get_token_value(from, pos, "energy_min", energy_min, energy_min);
        get_token_value(from, pos, "energy_max", energy_max, energy_max);
        get_token_value(from, pos, "rmsd_tol_uniq", rmsd_tol_uniq,
                        rmsd_tol_uniq);
        get_token_value(from, pos, "pop_size", pop_size, pop_size);
        get_token_value(from, pos, "max_iter", max_iter, max_iter);
        get_token_value(from, pos, "max_mut_tors", max_mut_tors, max_mut_tors);
        get_token_value(from, pos, "seed", seed, seed);
    }

    // Validate input:

    Assert::dynamic(dist_min >= 0.5, "bad dist_min < 0.5");
    Assert::dynamic(dist_max > dist_min, "bad dist_max <= dist_min");
    Assert::dynamic(rmsd_tol_uniq > 0.0, "bad rmsd_tol_uniq <= 0.0");
    Assert::dynamic(pop_size > 1, "bad pop_size <= 1");
    Assert::dynamic(max_iter > pop_size, "bad max_iter <= pop_size");
    Assert::dynamic(max_mut_tors >= 1, "bad max_mut_tors < 1");

    // Seed the random number engine:

    if (seed == 0) {
        std::random_device rd;
        std::seed_seq seed_seq_{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
        mt.seed(seed_seq_);
    }
    else {
        mt.seed(seed); // should only be used for testing purposes
    }

    // Initialize potential:

    pot.init(from);

    // Initialize population:

    init_population(to);
}

template <class Pot>
void Chem::Gamss<Pot>::solve(std::ostream& to)
{
    Stdutils::Format<char> line;
    line.width(52).fill('=');

    Stdutils::Format<double> fix;
    fix.fixed().width(12).precision(6);

    line.width(13).fill('-');
    to << "Local minima:\n" << line('-') << '\n';
    for (std::size_t i = 0; i < population.size(); ++i) {
        to << "Conformer: " << i + 1 << '\n'
           << "Energy: " << fix(population[i].energy) << '\n';
        Chem::Impl::print_geometry(to, mol.atoms(), population[i].xyz);
        to << '\n';
    }
}

//------------------------------------------------------------------------------

template <class Pot>
void Chem::Gamss<Pot>::init_population(std::ostream& to)
{
    Stdutils::Format<char> line;
    line.width(52).fill('=');

    int ipop = 0;
    int iter = 0;
    while (ipop < pop_size && iter < max_iter) {
        ++iter;

        // Generate new random structure with sensible geometry:
        Chem::Molecule m(mol);
        gen_rand_conformer(m);
        if (!geom_sensible(m)) {
            continue;
        }

        // Check if new random structure is blacklisted.
        if (is_blacklisted(m.cart_coord()) && iter > 0) {
            continue;
        }
        std::cout << m.cart_coord() << std::endl;

        // Add starting structure for local optimization to blacklist:
        blacklist.push_back(Chem::Conformer(m.elec_energy(), m.cart_coord()));

        // Perform local optimization:
        pot.run(m);

        if (m.elec_energy() >= energy_min && m.elec_energy() < energy_max) {
            // Add optimized structure to blacklist and population:
            blacklist.push_back(
                Chem::Conformer(m.elec_energy(), m.cart_coord()));
            population.push_back(
                Chem::Conformer(m.elec_energy(), m.cart_coord()));
            ++ipop;
        }
    }
    std::sort(population.begin(), population.end());

    to << "Genetic Algorithm Molecular Structure Search (GAMSS)\n"
       << line('=') << '\n';
    line.width(19).fill('-');
    to << "Initial population:\n" << line('-') << '\n';

    Stdutils::Format<double> fix;
    fix.fixed().width(12).precision(6);

    for (std::size_t i = 0; i < population.size(); ++i) {
        to << "Conformer: " << i + 1 << '\n'
           << "Energy: " << fix(population[i].energy) << '\n';
        Chem::Impl::print_geometry(to, mol.atoms(), population[i].xyz);
        to << '\n';
    }
}

template <class Pot>
void Chem::Gamss<Pot>::gen_rand_conformer(Chem::Molecule& m)
{
    int iter = 0;
    while (iter < max_mut_tors) {
        // Select a random dihedral angle:
        auto moiety = select_rand_dihedral(m);

        // Apply random variation to dihedral angle:
        std::uniform_real_distribution<> rnd_uni_real(-179.0, 180.0);
        double delta = rnd_uni_real(mt);
        m.int_coord().rotate_moiety(moiety, delta);

        ++iter;
    }
}

template <class Pot>
std::vector<int> Chem::Gamss<Pot>::select_rand_dihedral(const Chem::Molecule& m)
{
    std::vector<Numlib::Vec<int>> connect = m.int_coord().get_connectivities();
    std::uniform_int_distribution<> rnd_uni_int(2, connect.size() - 1);
    int index = rnd_uni_int(mt);
    auto dihedral = connect[index];

    std::vector<int> res;
    for (std::size_t i = 2; i < connect.size(); ++i) {
        if (connect[i] == dihedral) {
            res.push_back(i);
        }
    }
    return res;
}

template <class Pot>
bool Chem::Gamss<Pot>::geom_sensible(const Chem::Molecule& m) const
{
    bool geom_ok = true;
    Numlib::Mat<double> dist_mat;
    Numlib::pdist_matrix(dist_mat, m.cart_coord());
    for (auto v : dist_mat) {
        if (v > 0.0 && v < dist_min) { // avoid too close atoms
            geom_ok = false;
            break;
        }
    }
    if (geom_ok) { // avoid too long bond distances
        for (std::size_t i = 0; i < m.num_atoms(); ++i) {
            if (m.int_coord().get_distance(i) >= dist_max) {
                geom_ok = false;
                break;
            }
        }
    }
    return geom_ok;
}

template <class Pot>
bool Chem::Gamss<Pot>::is_blacklisted(const Numlib::Mat<double>& xyz) const
{
    bool res = false;
    for (const auto& conformer : blacklist) {
        if (Numlib::kabsch_rmsd(conformer.xyz, xyz) < rmsd_tol_uniq) {
            res = true;
            break;
        }
    }
    return res;
}

template class Chem::Gamss<Chem::Gaussian>;
template class Chem::Gamss<Chem::Mopac>;

