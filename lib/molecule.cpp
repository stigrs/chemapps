// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/molecule.h>
#include <chem/io.h>
#include <stdutils/stdutils.h>

Chem::Molecule::Molecule(std::istream& from,
                         std::ostream& to,
                         const std::string& key,
                         bool verbose)
    : elec_(from, key),
      geom_(from, key),
      rot_(from, key, geom_.atoms(), geom_.get_xyz()),
      vib_(from,
           key,
           geom_.atoms(),
           rot_.get_xyz_paxis(),
           rot_.principal_axes()),
      tor_(from,
           key,
           geom_.atoms(),
           rot_.get_xyz_paxis(),
           rot_.principal_axes(),
           rot_.principal_moments())
{
    if (verbose) {
        Stdutils::Format<char> line;
        line.width(15 + key.size()).fill('=');

        Stdutils::Format<double> fix;
        fix.fixed().precision(6);

        to << "Input data on " << key << ":\n" << line('=') << '\n';
        to << "Electronic energy: " << fix(elec_.energy()) << " Hartree\n"
           << "Charge: " << elec_.charge() << '\n'
           << "Spin multiplicity: " << elec_.spin_mult() << '\n';
        Chem::print_spin_orbit_states(to, elec_.spin_orbit_degen(),
                                      elec_.spin_orbit_energy());
        to << "\nInput orientation:\n";
        Chem::print_geometry(to, geom_.atoms(), geom_.get_xyz());
        Chem::print_atomic_masses(to, geom_.atoms());
        vib_.print(to);
    }
}

double Chem::Molecule::tot_mass() const
{
    double res = 0.0;
    for (const auto& at : atoms()) {
        res += at.atomic_mass;
    }
    return res;
}

