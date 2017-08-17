/**
   @file molecule.cpp

   This file is part of ChemApps - A C++ Chemistry Toolkit

   Copyright (C) 2016-2017  Stig Rune Sellevag

   ChemApps is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   ChemApps is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <chem/input.h>
#include <chem/molecule.h>
#include <chem/molecule_io.h>
#include <chem/utils.h>
#include <map>

Molecule::Molecule(const Molecule& mol)
    : atoms(mol.atoms), xyz(mol.xyz), elec_state(mol.elec_state)
{
    title = mol.title;
    geom_unit = mol.geom_unit;
    charge = mol.charge;
    elec_energy = mol.elec_energy;

    zmat = std::make_shared<Zmatrix>(*mol.zmat);
    rot = std::make_shared<Molrot>(*mol.rot);
    vib = std::make_shared<Molvib>(*mol.vib);
    tor = std::make_shared<Torsion>(*mol.tor);
}

void Molecule::print_data(std::ostream& to, const std::string& key) const
{
    chem::Format<char> line;
    line.width(15 + key.size()).fill('=');

    chem::Format<double> fix;
    fix.fixed().precision(6);

    to << "Input data on " << key << ":\n" << line('=') << '\n';

    chem::print_elec_states(to, elec_state);
    to << "Electronic energy: " << fix(elec_energy) << " Hartree\n"
       << "Charge: " << charge << "\n\n"
       << "Input orientation:\n";
    chem::print_geometry(to, atoms, xyz, geom_unit);
    chem::print_atomic_masses(to, atoms);
    vib->print(to);
}

void Molecule::init(std::istream& from,
                    std::ostream& to,
                    const std::string& key,
                    bool verbose)
{
    // Read input data:

    arma::vec elec_state_def(2);
    elec_state_def(0) = 1.0;
    elec_state_def(1) = 0.0;

    std::map<std::string, Input> input_data;
    input_data["geom_unit"]   = Input(geom_unit, "angstrom");
    input_data["charge"]      = Input(charge, 0);
    input_data["elec_state"]  = Input(elec_state, elec_state_def);
    input_data["elec_energy"] = Input(elec_energy, 0.0);

    bool found = chem::find_section(from, key);
    if (found) {
        std::string token;
        while (from >> token) {
            if (token == "End") {
                break;
            }
            else if (token == "geometry") {
                chem::read_xyz_format(from, atoms, xyz, title);
            }
            else {
                auto it = input_data.find(token);
                if (it != input_data.end()) {
                    from >> it->second;
                }
            }
        }
    }
    else {
        throw Mol_error("cannot find " + key + " section");
    }

    // Check if initialized:

    for (auto it = input_data.begin(); it != input_data.end(); ++it) {
        if (!it->second.is_init()) {
            throw Mol_error(it->first + " not initialized");
        }
    }

    // Initialize Z matrix object:

    zmat = std::make_shared<Zmatrix>(atoms, xyz);

    // Initialize molecular rotations object:

    rot = std::make_shared<Molrot>(from, key, atoms, xyz);

    // Initialize molecular vibrations object:

    vib = std::make_shared<Molvib>(from, key);

    // Initialize molecular torsions object:

    tor = std::make_shared<Torsion>(from, key, *rot);

    // Write input data to output stream:

    if (verbose) {
        print_data(to, key);
    }
}
