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
    : geom_(from, key)
{
    rot_ = Chem::Rotation(from, key, geom_.atoms(), geom_.get_xyz());
    vib_ = Chem::Vibration(from, key, geom_.atoms(), rot_.get_xyz_paxis(),
                           rot_.principal_axes());

    if (verbose) {
        Stdutils::Format<char> line;
        line.width(15 + key.size()).fill('=');

        Stdutils::Format<double> fix;
        fix.fixed().precision(6);

        to << "Input data on " << key << ":\n" << line('=') << '\n';
#if 0
        to << "Electronic energy: " << fix(elec.elec_energy()) << " Hartree\n"
           << "Charge: " << elec.net_charge() << '\n'
           << "Spin multiplicity: " << elec.spin_mult() << '\n';
        Chem::Impl::print_spin_orbit_states(to, elec.spin_orbit_degen(),
                                            elec.spin_orbit_energy());
#endif
        to << "\nInput orientation:\n";
        Chem::print_geometry(to, geom_.atoms(), geom_.get_xyz());
        Chem::print_atomic_masses(to, geom_.atoms());
        vib_.print(to);
    }
}

