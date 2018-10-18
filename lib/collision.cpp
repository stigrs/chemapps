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

#include <chem/collision.h>
#include <chem/impl/molecule_io.h>
#include <chem/periodic_table.h>
#include <stdutils/stdutils.h>
#include <stdexcept>
#include <cassert>

Chem::Collision::Collision(std::istream& from, const std::string& key)
{
    using namespace Stdutils;

    // Populate local sigma and epsilon values:
    set_sigma_local_values();
    set_epsilon_local_values();

    // Read input data:
    std::string coll_int_str;

    auto pos = find_token(from, key);
    if (pos != -1) {
        get_token_value(from, pos, "coll_integral", coll_int_str, "forst");
        get_token_value(from, pos, "mass_bath", mass_bath);
        get_token_value(from, pos, "mass_mol", mass_mol);
        get_token_value(from, pos, "epsilon_bath", epsilon_bath);
        get_token_value(from, pos, "epsilon_mol", epsilon_mol);
        get_token_value(from, pos, "sigma_bath", sigma_bath);
        get_token_value(from, pos, "sigma_mol", sigma_mol);
        get_token_value(from, pos, "vibr_high", vibr_high, 0.0);
        pos = find_token(from, "mol_formula");
        if (pos != -1) {
            Chem::Impl::read_mol_formula(from, mol_formula);
        }
    }
    // Validate input data:

    if (coll_int_str == "troe") {
        coll_integral = troe;
    }
    else if (coll_int_str == "forst") {
        coll_integral = forst;
    }
    else {
        throw std::runtime_error("bad coll_integral: " + coll_integral_str);
    }
    assert(mass_bath > 0.0);
    assert(mass_mol > 0.0);
    assert(epsilon_bath > 0.0);
    assert(epsilon_mol > 0.0);
    assert(sigma_bath > 0.0);
    assert(sigma_mol > 0.0);
}

void Collision::biased_random_walk(double temp, std::ostream& to) const
{
    Stdutils::Format<char> line;
    line.width(35).fill('-');

    to << "Lennard-Jones collision parameters:\n" << line('-') << '\n';
    to << "Collision integral:\t\t";
    switch (coll_integral) {
    case troe:
        to << "Troe\n";
        break;
    case forst:
        to << "Forst\n";
        break;
    }

    to << "Collision diameter bath gas:\t" << sigma_bath << " angstrom\n"
       << "Collision diameter molecule:\t" << sigma_mol << " angstrom\n"
       << "Collision diameter complex:\t" << sigma_complex() << " angstrom\n"
       << "Collision well depth bath gas:\t" << epsilon_bath << " K\n"
       << "Collision well depth molecule:\t" << epsilon_mol << " K\n"
       << "Collision well depth complex:\t" << epsilon_complex() << " K\n"
       << "Mass bath gas:\t\t\t" << mass_bath << " amu\n"
       << "Mass molecule:\t\t\t" << mass_mol << " amu\n"
       << "Reduced mass complex:\t\t" << reduced_mass() << " amu\n";

    to << "Highest vibrational frequency:\t" << vibr_high << " cm-1\n"
       << "Average atom/atom mass:\t\t" << average_mass() << " amu\n"
       << "Local well depth complex:\t" << epsilon_local() << " K\n"
       << "Local collision diam. complex:\t" << sigma_local() << " angstrom\n";

    line.width(27).fill('-');
    to << "\nBiased random walk results:\n" << line('-') << '\n';

    double om22 = coll_omega22(temp);
    double d = dist_interact(temp);
    double b = impact_parameter(temp);
    double etr = energy_trans_avg(temp);
    double tc = collision_time(temp);
    double a = a_decay_parameter(temp);
    double c = c_autocorr_osc_freq();
    double edot = mean_sqr_int_energy_change(temp);
    double s = s_parameter(temp);
    double e2 = std::sqrt(mean_sqr_energy_transfer_coll(temp));

    to << "Collision integral:\t\t" << om22 << '\n'
       << "Closest interaction distance:\t" << d << " angstrom\n"
       << "Impact parameter:\t\t" << b << " angstrom\n"
       << "Average translational energy:\t" << etr << " cm-1\n"
       << "Collision time:\t\t\t" << tc << " s\n"
       << "A decay parameter:\t\t" << a << " s-1\n"
       << "Autocorr. osc. frequency (C):\t" << c << " s-1\n"
       << "<Edot(i)^2>:\t\t\t" << edot << " cm-2 s-2\n"
       << "BRW parameter (s):\t\t" << s << " cm-1\n"
       << "sqrt(<E^2>):\t\t\t" << e2 << " cm-1\n";
}

double Collision::average_mass() const
{
    double mavg = 0.0; // average mass of atoms in molecule

    int natoms = 0;
    if (!mol_formula.empty()) {
        for (auto& i : mol_formula) {
            natoms += i.stoich;
            mavg += i.stoich * Periodic_table::get_atomic_mass(i.atom);
        }
        mavg /= natoms;
    }
    // Eq. 35b in Lim and Gilbert (1990)
    return 1.0 / ((1.0 / (mavg * natoms - mavg)) + (1.0 / mavg));
}

double Collision::sigma_local() const
{
    using namespace Periodic_table;

    double sloc = 0.0;

    int natoms = 0;
    for (auto& i : mol_formula) {
        natoms += i.stoich;
        sloc += i.stoich * sigma_loc_val[get_atomic_number(i.atom)];
    }
    sloc /= natoms;
    return 0.5 * (sloc + sigma_bath);
}

double Collision::epsilon_local() const
{
    using namespace Periodic_table;

    double eloc = 0.0;

    int natoms = 0;
    for (auto& i : mol_formula) {
        natoms += i.stoich;
        eloc += i.stoich * epsilon_loc_val[get_atomic_number(i.atom)];
    }
    eloc /= natoms;
    return std::sqrt(eloc * epsilon_bath);
}

double Collision::coll_omega22(double temp) const
{
    double red_temp = temp / epsilon_complex();

    double omega22 = 0.0;
    switch (coll_integral) {
    case forst: // eq. A4.10 in Forst (2003)
        omega22 = 1.16145 / std::pow(red_temp, 0.148774) +
                  0.52487 / std::exp(0.7732 * red_temp) +
                  2.16178 / std::exp(2.437887 * red_temp);
        break;
    case troe: // eq. 3.2 in Troe (1977)
        omega22 = 1.0 / (0.697 + 0.5185 * std::log10(red_temp));
        break;
    default: // unknown
        omega22 = 0.0;
    }
    return omega22;
}

double Collision::collision_time(double temp) const
{
    const double ang2m = 1.0e-10; // angstrom to meters
    const double dr = 0.0005;     // step size

    const double d = dist_interact(temp);
    const double b = impact_parameter(temp);
    const double etr = energy_trans_avg(temp) * datum::icm_to_K;
    const double eps = epsilon_complex();
    const double sig = sigma_complex();

    double r = d;
    double veff = 0.0;
    double tc = 0.0;

    while (veff < etr) { // integrate eq. 32 in Lim and Gilbert (1990)
        double sigr6 = std::pow(sig / r, 6.0);
        double sigr12 = sigr6 * sigr6;
        veff = 4.0 * eps * (sigr12 - sigr6) + etr * std::pow(b / r, 2.0);
        if (veff >= etr) {
            break;
        }
        else {
            tc += 1.0 / std::sqrt(etr - veff);
            r -= dr;
        }
    }
    return tc *= (std::sqrt(2.0 * reduced_mass() * datum::m_u) * dr * ang2m /
                  std::sqrt(datum::k));
}

void Collision::set_sigma_local_values()
{
    using namespace Periodic_table;

    sigma_loc_val.resize(get_max_atomic_number(), 0.0);
    sigma_loc_val[get_atomic_number("H")] = 3.0;
    sigma_loc_val[get_atomic_number("C")] = 3.2;
    sigma_loc_val[get_atomic_number("N")] = 3.2;
    sigma_loc_val[get_atomic_number("O")] = 3.2;
    sigma_loc_val[get_atomic_number("S")] = 3.4;
    sigma_loc_val[get_atomic_number("F")] = 3.2;
    sigma_loc_val[get_atomic_number("Cl")] = 3.4;
    sigma_loc_val[get_atomic_number("Br")] = 3.6;
    sigma_loc_val[get_atomic_number("I")] = 4.0;
}

void Collision::set_epsilon_local_values()
{
    using namespace Periodic_table;

    epsilon_loc_val.resize(get_max_atomic_number(), 0.0);
    epsilon_loc_val[get_atomic_number("H")] = 6.5;
    epsilon_loc_val[get_atomic_number("C")] = 20.3;
    epsilon_loc_val[get_atomic_number("N")] = 20.3;
    epsilon_loc_val[get_atomic_number("O")] = 20.3;
    epsilon_loc_val[get_atomic_number("S")] = 120.0;
    epsilon_loc_val[get_atomic_number("F")] = 20.3;
    epsilon_loc_val[get_atomic_number("Cl")] = 120.0;
    epsilon_loc_val[get_atomic_number("Br")] = 190.0;
    epsilon_loc_val[get_atomic_number("I")] = 230.0;
}

double Collision::a_decay_parameter(double temp) const
{
    using namespace Numlib::Constants;

    const double ang2m = 1.0e-10;    // angstrom to meters
    const double dr = 0.005 * ang2m; // step size

    const double ebar = energy_trans_avg(temp) * icm_to_K * k;
    const double eps = epsilon_local() * k;
    const double sig = sigma_local() * ang2m;
    const double d = std::max(sig, sig * coll_omega22(temp));
    const double b = (2.0 / 3.0) * d;

    double r = d;
    double sigr6 = 1.0;
    double sigr12 = 1.0;
    double veff = 0.0;
    double a = 0.0;

    while (veff < ebar) { // evaluate eq. 36 in Lim and Gilbert (1990)
        sigr6 = std::pow(sig / r, 6.0);
        sigr12 = sigr6 * sigr6;
        veff = 4.0 * eps * (sigr12 - sigr6) + ebar * std::pow(b / r, 2.0);
        if (veff >= ebar) {
            break;
        }
        else {
            r -= dr;
        }
    }
    // Evaluate eq. 34b in Lim and Gilbert (1990):

    a = std::abs((4.0 * eps * (-12.0 * sigr12 + 6.0 * sigr6) -
                  2.0 * ebar * std::pow(b / r, 2.0)) /
                 r);
    return a /= std::sqrt(0.5 * ebar * average_mass() * m_u);
}

double Collision::mean_sqr_int_energy_change(double temp) const
{
    // Evaluate eq. 38 in Lim and Gilbert (1990):

    const double ang2m = 1.0e-10; // angstrom to meter

    const double ebar = energy_trans_avg(temp) * datum::icm_to_K * datum::k;
    const double eps = epsilon_local() * datum::k;
    const double sig = sigma_local() * ang2m;
    const double mlight = mol_mass_lightest() * datum::m_u;
    const double nu = vibr_high * datum::c_0 * 100.0;
    const double k = 4.0 * datum::pi * datum::pi * mlight * nu * nu;
    const double deltax = std::sqrt(2.0 * ebar / k);
    const double f = 6.0 / 5.0;
    const double x = sig / f;

    double deltav = 4.0 * eps * (deltax / x) *
                    (-12.0 * std::pow(f, 12.0) + 6.0 * std::pow(f, 6.0));

    // Evaluate eq. 40 in Lim and Gilbert (1990):

    return std::pow((ebar - std::min(std::abs(deltav), 0.5 * ebar)) * nu, 2.0) /
           (datum::J_to_icm * datum::J_to_icm);
}

double Collision::mol_mass_lightest() const
{
    double mlight = 1.0e+6;
    for (auto& i : mol_formula) {
        double mtmp = ptable::get_atomic_mass(i.atom);
        if (mtmp < mlight) {
            mlight = mtmp;
        }
    }
    return mlight;
}
