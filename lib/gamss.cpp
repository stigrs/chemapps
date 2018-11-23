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

template <class Pot>
Chem::Gamss<Pot>::Gamss(std::istream& from) : mol(from)
{
    using namespace Stdutils;
    const std::string key = "Gamss";

    // Read input:

    dist_min = 0.5;
    pop_size = 100;
    max_trials = 20;
    seed = 0;

    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "dist_min", dist_min);
        get_token_value(from, pos, "pop_size", pop_size);
        get_token_value(from, pos, "max_trials", max_trials);
        get_token_value(from, pos, "seed", seed);
    }

    // Validate input:

    Assert::dynamic(dist_min > 0.1, "bad dist_min <= 0.1");
    Assert::dynamic(pop_size > 1, "bad pop_size <= 1");
    Assert::dynamic(max_trials > 1, "bad max_trials <= 1");

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

    init_population();
}

template <class Pot>
void Chem::Gamss<Pot>::solve(std::ostream& to)
{
    Stdutils::Format<char> line;
    line.width(44).fill('=');

    Stdutils::Format<double> fix;
    fix.fixed().width(12).precision(6);

    to << "Genetic Algorithm Molecular Structure Search\n" << line('=') << '\n';

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
void Chem::Gamss<Pot>::init_population()
{
    for (int i = 0; i < pop_size; ++i) {
        // Generate new random structure with sensible geometry:
        Chem::Molecule m(mol);
        gen_rand_conformer(m);

        // Add starting structure for local optimization to blacklist:
        blacklist.push_back(Chem::Conformer(m.elec_energy(), m.cart_coord()));

        // Perform local optimization:
        pot.run(m);

        // Add optimized structure to blacklist and population:
        blacklist.push_back(Chem::Conformer(m.elec_energy(), m.cart_coord()));
        population.push_back(Chem::Conformer(m.elec_energy(), m.cart_coord()));
    }
}

template <class Pot>
void Chem::Gamss<Pot>::gen_rand_conformer(Molecule& m)
{
    int iter = 0;
    while (iter < max_trials) {
        // Select a random dihedral angle:
        auto moiety = select_rand_dihedral(m);

        // Apply random variation to dihedral angle:
        std::uniform_real_distribution<> rnd_uni_real(-180.0, 180.0);
        double delta = rnd_uni_real(mt);
        m.int_coord().rotate_moiety(moiety, delta);

        // Check if geometry is sensible:
        if (geom_sensible(m)) {
            break;
        }
        ++iter;
    }
    if (iter >= max_trials) {
        throw std::runtime_error(
            "could not generate random conformer, increase max_trials");
    }
}

template <class Pot>
std::vector<int> Chem::Gamss<Pot>::select_rand_dihedral(const Chem::Molecule& m)
{
    std::vector<Numlib::Vec<int>> connect = m.int_coord().get_connectivities();
    for (auto v : connect) {
        std::cout << v << std::endl;
    }
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
        }
    }
    return geom_ok;
}

template class Chem::Gamss<Chem::Gaussian>;
template class Chem::Gamss<Chem::Mopac>;

