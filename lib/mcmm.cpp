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

#include <chem/gaussian.h>
#include <chem/mcmm.h>
#include <chem/molecule_io.h>
#include <chem/mopac.h>
#include <srs/utils.h>
#include <limits>
#include <map>
#include <numeric>


template <class Pot>
Mcmm<Pot>::Mcmm(std::istream& from,
                Molecule& mol_,
                const std::string& key,
                bool verbose_)
    : mol(mol_), verbose(verbose_)
{
    // Read input data:

    const double emin_def = -std::numeric_limits<double>::max();

    std::map<std::string, srs::Input> input_data;
    input_data["xtol"]      = srs::Input(xtol, 5.0e-3);
    input_data["etol"]      = srs::Input(etol, 1.0e-4);
    input_data["emin"]      = srs::Input(emin, emin_def);
    input_data["emax"]      = srs::Input(emax, 0.0);
    input_data["rmin"]      = srs::Input(rmin, 0.7414);  // experimental r(H-H)
    input_data["temp"]      = srs::Input(temp, 298.15);
    input_data["maxiter"]   = srs::Input(maxiter, 1000);
    input_data["miniter"]   = srs::Input(miniter, 100);
    input_data["maxreject"] = srs::Input(maxreject, 100);
    input_data["nminima"]   = srs::Input(nminima, 20);
    input_data["seed"]      = srs::Input(seed, 0);

    bool found = srs::find_section(from, key);
    if (found) {
        std::string token;
        while (from >> token) {
            if (token == "End") {
                break;
            }
            else {
                auto it = input_data.find(token);
                if (it != input_data.end()) {
                    from >> it->second;
                }
            }
        }
    }

    // Check if initialized:

    for (auto it = input_data.begin(); it != input_data.end(); ++it) {
        if (!it->second.is_init()) {
            throw Mcmm_error(it->first + " not initialized");
        }
    }

    // Validate input:

    if (xtol <= 0.0) {
        throw Mcmm_error("bad xtol <= 0.0");
    }
    if (etol <= 0.0) {
        throw Mcmm_error("bad etol <= 0.0");
    }
    if (rmin <= 0.0) {
        throw Mcmm_error("bad rmin <= 0.0");
    }
    if (temp <= 0.0) {
        throw Mcmm_error("bad temp <= 0.0");
    }
    if (maxiter < nminima) {
        throw Mcmm_error("bad maxiter < nminima");
    }
    if (maxiter < miniter) {
        throw Mcmm_error("bad maxiter < miniter");
    }
    if (maxreject < 1) {
        throw Mcmm_error("bad maxreject < 1");
    }
    if (nminima < 1) {
        throw Mcmm_error("bad nminima < 1");
    }

    // Initialize potential:

    pot.init(from);

    // Initialize iterators:

    kiter   = 0;
    nreject = 0;
    naccept = 0;

    // Initialize storage containers:

    xcurr = mol.get_xyz();
    ecurr = mol.get_elec_energy();

    // Seed the random number engine:

    if (seed == 0) {
        std::random_device rd;
        std::seed_seq seed_seq_{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
        mt.seed(seed_seq_);
    }
    else {
        mt.seed(seed);  // should only be used for testing purposes
    }
}

template <class Pot>
void Mcmm<Pot>::solve(std::ostream& to)
{
    global_min_found = false;
    while (!global_min_found) {
        new_conformer();
        update();
        if (check_exit()) {
            sort_conformers();
            global_min_found = true;
        }
    }
    if (verbose) {
        double eglobal_min = *std::min_element(eglobal.begin(), eglobal.end());
        srs::Format<char> line;
        line.width(41).fill('=');
        srs::Format<int> ifix;
        ifix.fixed().width(5);
        srs::Format<double> dfix;
        dfix.fixed().width(8).precision(2);
        to << "Monte Carlo Multiple Minima (MCMM) Solver\n"
           << line('=') << '\n';
        to << "Temperature:\t" << dfix(temp) << '\n'
           << "Iterations:\t" << ifix(kiter) << " out of " << maxiter << '\n'
           << "Rejections:\t" << ifix(nreject) << " out of " << maxreject
           << "\n\n";
        dfix.fixed().width(12).precision(6);
        line.width(15).fill('-');
        to << "Global minimum:\n"
           << line('-') << '\n'
           << "Energy: " << dfix(eglobal_min) << '\n';
        chem::print_geometry(to, mol.get_atoms(), xglobal);
        to << '\n';

        line.width(13).fill('-');
        to << "Local minima:\n" << line('-') << '\n';
        for (std::size_t i = 0; i < conformers.size(); ++i) {
            to << "Conformer: " << i + 1 << '\n'
               << "Energy: " << dfix(conformers[i].energy) << '\n';
            chem::print_geometry(to, mol.get_atoms(), conformers[i].xyz);
            to << '\n';
        }
    }
}

template <class Pot>
bool Mcmm<Pot>::check_exit() const
{
    bool finished = false;
    if (ecurr < emin) {
        finished = true;
    }
    if (kiter >= maxiter) {
        finished = true;
    }
    if (nreject >= maxreject) {
        finished = true;
    }
    std::vector<double> ediff(eglobal.size());
    std::adjacent_difference(eglobal.begin(), eglobal.end(), ediff.begin());
    if (ediff.size() > 1) {
        double ediff_max = *std::max_element(ediff.begin(), ediff.end());
        if ((ediff_max < etol) && (kiter >= miniter)) {
            finished = true;
        }
    }
    return finished;
}

template <class Pot>
bool Mcmm<Pot>::accept_energy(double enew)
{
    bool accept  = false;
    double ediff = enew - ecurr;
    if (enew > emax) {
        nreject += 1;
    }
    else if (ediff < 0.0) {
        accept = true;
    }
    else {
        double h = std::exp(-ediff / temp);
        std::uniform_real_distribution<> rnd_real_uni(0.0, 1.0);
        if (h > rnd_real_uni(mt)) {
            accept = true;
        }
        else {
            nreject += 1;
        }
    }
    return accept;
}

template <class Pot>
bool Mcmm<Pot>::accept_geom_dist(const Molecule& m) const
{
    bool geom_ok = true;
    srs::dmatrix dist_mat;
    srs::pdist_matrix(dist_mat, m.get_xyz());
    for (auto v : dist_mat) {
        if (v > 0.0 && v < rmin) {  // avoid too close atoms
            geom_ok = false;
        }
    }
    return geom_ok;
}

template <class Pot>
bool Mcmm<Pot>::duplicate(const Molecule& m) const
{
    bool duplicate = false;
    for (std::size_t i = 0; i < conformers.size(); ++i) {  // check geometry
        duplicate = srs::rmsd(conformers[i].xyz, m.get_xyz()) <= xtol;
        if (duplicate) {  // check energy
            double ediff = std::abs(conformers[i].energy - m.get_elec_energy());
            duplicate    = ediff < etol;
        }
    }
    return duplicate;
}

template <class Pot>
void Mcmm<Pot>::new_conformer()
{
    const unsigned ntrials = 20;
    // Generate a new random conformer by using the uniform usage scheme:
    for (unsigned i = 0; i < ntrials; ++i) {
        srs::dmatrix xnew = mol.get_xyz();
        uniform_usage(xnew);
        mol.set_xyz(xnew);
        gen_rand_conformer(mol);
        if (accept_geom_dist(mol)) {  // check geometry constraints
            break;
        }
    }

    // Perform geometry optimization:
    pot.run(mol);

    // Check acceptance:
    if (accept_energy(mol.get_elec_energy())) {
        if (!duplicate(mol)) {  // store new conformer
            xcurr = mol.get_xyz();
            ecurr = mol.get_elec_energy();
            save_conformer(mol);
            naccept += 1;
        }
    }
}

template <class Pot>
void Mcmm<Pot>::update()
{
    kiter += 1;

    // Update global minimum if appropriate:
    if (!eglobal.empty()) {
        if (ecurr <= *std::min_element(eglobal.begin(), eglobal.end())) {
            eglobal.push_back(ecurr);
            xglobal = xcurr;
        }
    }
    else {
        eglobal.push_back(ecurr);
        xglobal = xcurr;
    }
}

template <class Pot>
void Mcmm<Pot>::uniform_usage(srs::dmatrix& xnew)
{
    xnew = xcurr;
    if (!conformers.empty()) {
        std::vector<Conformer> index;
        int istart     = 0;
        int min_nstart = kiter;
        // Find least used:
        for (std::size_t i = 0; i < conformers.size(); ++i) {
            if (conformers[i].iter <= min_nstart) {
                index.push_back(conformers[i]);
            }
        }
        if (!conformers.empty()) {  // find least used with lowest energy
            auto res = std::max_element(conformers.begin(), conformers.end());
            auto it  = std::distance(conformers.begin(), res);
            double emin_ = conformers[it].energy;
            for (std::size_t i = 0; i < index.size(); ++i) {
                if (index[i].energy < emin_) {
                    emin_  = index[i].energy;
                    istart = i;
                }
            }
        }
        conformers[istart].iter += 1;
        xnew = conformers[istart].xyz;
    }
}

template <class Pot>
inline void Mcmm<Pot>::sort_conformers()
{
    std::sort(conformers.begin(), conformers.end());
    if (conformers.size() > nminima) {
        std::vector<Conformer> tmp;
        for (std::size_t i = 0; i < nminima; ++i) {
            tmp.push_back(Conformer(conformers[i]));
        }
        conformers.resize(nminima);
        conformers = tmp;
    }
}

template <class Pot>
std::vector<int> Mcmm<Pot>::select_rand_dihedral(const Molecule& m)
{
    std::vector<srs::ivector> connect = m.get_zmat().get_connectivities();
    std::uniform_int_distribution<> rnd_uni_int(2, connect.size() - 1);
    int index             = rnd_uni_int(mt);
    srs::ivector dihedral = connect[index];

    std::vector<int> res(0);
    for (std::size_t i = 2; i < connect.size(); ++i) {
        if (connect[i] == dihedral) {
            res.push_back(i);
        }
    }
    return res;
}

template class Mcmm<Gaussian>;
template class Mcmm<Mopac>;
