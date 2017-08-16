/**
   @file molecule_io.cpp

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

#include <chem/molecule_io.h>
#include <chem/ptable.h>
#include <chem/utils.h>
#include <sstream>

void chem::read_xyz_format(std::istream& from,
                           std::vector<Element>& atoms,
                           arma::mat& xyz,
                           std::string& title)
{
    // Get number of atoms:
    int natoms;
    from >> natoms;
    from.ignore();  // need to consume '\n' before reading title line
    from.clear();

    xyz.set_size(natoms, 3);
    atoms.resize(natoms);

    std::string symbol;
    double x;
    double y;
    double z;

    // Read title line:
    std::getline(from, title);
    title = chem::trim(title, " ");

    /// Read XYZ coordinates:
    for (int i = 0; i < natoms; ++i) {
        from >> symbol >> x >> y >> z;
        atoms[i] = ptable::get_element(symbol);
        xyz(i, 0) = x;
        xyz(i, 1) = y;
        xyz(i, 2) = z;
    }
}

void chem::read_zmat_format(std::istream& from,
                            std::vector<Element>& atoms,
                            arma::vec& distances,
                            arma::vec& angles,
                            arma::vec& dihedrals,
                            arma::ivec& bond_connect,
                            arma::ivec& angle_connect,
                            arma::ivec& dihedral_connect)
{
    distances.clear();
    angles.clear();
    dihedrals.clear();
    bond_connect.clear();
    angle_connect.clear();
    dihedral_connect.clear();
    atoms.clear();

    from.clear();
    from.seekg(0, std::ios_base::beg);

    std::string line;
    std::string symbol;
    int natoms = 0;
    int pos    = 0;

    // Count number of atoms:

    while (std::getline(from, line)) {
        std::istringstream iss(line);
        iss >> symbol;
        if (symbol.find("#") != std::string::npos) {
            continue;
        }
        else if (symbol.find("zmatrix") != std::string::npos) {
            pos = from.tellg();
        }
        else if (ptable::atomic_symbol_is_valid(symbol)) {
            atoms.push_back(ptable::get_element(symbol));
            natoms += 1;
        }
        else if (symbol.empty()) {
            break;
        }
    }

    // Load Z matrix:

    from.clear();
    from.seekg(pos, std::ios_base::beg);

    distances        = arma::zeros<arma::vec>(natoms);
    angles           = arma::zeros<arma::vec>(natoms);
    dihedrals        = arma::zeros<arma::vec>(natoms);
    bond_connect     = arma::zeros<arma::ivec>(natoms);
    angle_connect    = arma::zeros<arma::ivec>(natoms);
    dihedral_connect = arma::zeros<arma::ivec>(natoms);

    if (natoms > 0) {
        std::getline(from, line);  // first atom is already read
    }
    int iat1;
    double distance;
    if (natoms > 1) {
        std::getline(from, line);
        std::istringstream iss(line);
        iss >> symbol >> iat1 >> distance;
        distances[1]    = distance;
        bond_connect[1] = iat1 - 1;
    }
    int iat2;
    double angle;
    if (natoms > 2) {
        std::getline(from, line);
        std::istringstream iss(line);
        iss >> symbol >> iat1 >> distance >> iat2 >> angle;
        distances[2]     = distance;
        bond_connect[2]  = iat1 - 1;
        angles[2]        = angle;
        angle_connect[2] = iat2 - 1;
    }
    int iat3;
    double dihedral;
    if (natoms > 3) {
        for (int i = 3; i < natoms; ++i) {
            std::getline(from, line);
            std::istringstream iss(line);
            // clang-format off
            iss >> symbol   >> iat1 
                >> distance >> iat2 
                >> angle    >> iat3 
                >> dihedral;
            // clang-format on
            distances[i]        = distance;
            bond_connect[i]     = iat1 - 1;
            angles[i]           = angle;
            angle_connect[i]    = iat2 - 1;
            dihedrals[i]        = dihedral;
            dihedral_connect[i] = iat3 - 1;
        }
    }
}

void chem::print_xyz_format(std::ostream& to,
                            const std::vector<Element>& atoms,
                            const arma::mat& xyz,
                            const std::string& title)
{
    chem::Format<double> fix;
    fix.fixed().width(10);

    to << atoms.size() << '\n';
    to << title << '\n';
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        to << atoms[i].atomic_symbol << '\t';
        for (arma::uword j = 0; j < xyz.n_cols; ++j) {
            to << fix(xyz(i, j)) << '\t';
        }
        to << '\n';
    }
}

void chem::print_zmat_format(std::ostream& to,
                             std::vector<Element>& atoms,
                             arma::vec& distances,
                             arma::vec& angles,
                             arma::vec& dihedrals,
                             arma::ivec& bond_connect,
                             arma::ivec& angle_connect,
                             arma::ivec& dihedral_connect)
{
    chem::Format<double> dfix;
    dfix.fixed().width(10).precision(4);

    chem::Format<int> ifix;
    ifix.fixed().width(3);

    if (!atoms.empty()) {
        to << atoms[0].atomic_symbol << '\n';
    }
    if (atoms.size() > 1) {
        to << atoms[1].atomic_symbol << '\t' << ifix(bond_connect(1) + 1)
           << "  " << dfix(distances(1)) << '\n';
    }
    if (atoms.size() > 2) {
        to << atoms[2].atomic_symbol << '\t' << ifix(bond_connect(2) + 1)
           << "  " << dfix(distances(2)) << "  " << ifix(angle_connect(2) + 1)
           << "  " << dfix(angles(2)) << '\n';
    }
    if (atoms.size() > 3) {
        for (std::size_t i = 3; i < atoms.size(); ++i) {
            to << atoms[i].atomic_symbol << '\t' << ifix(bond_connect(i) + 1)
               << "  " << dfix(distances(i)) << "  "
               << ifix(angle_connect(i) + 1) << "  " << dfix(angles(i)) << "  "
               << ifix(dihedral_connect(i) + 1) << "  " << dfix(dihedrals(i))
               << '\n';
        }
    }
}

void chem::print_elec_states(std::ostream& to, const arma::vec& elec_state)
{
    chem::Format<char> line;
    line.width(34).fill('-');

    chem::Format<double> fix;
    fix.fixed().width(6).precision(2);

    to << "Electronic states:\n"
       << line('-') << '\n'
       << " #\tEnergy/cm^-1\tDegeneracy\n"
       << line('-') << '\n';
    int it = 1;
    for (arma::uword i = 0; i < elec_state.size(); i += 2) {
        to << " " << it << '\t' << fix(elec_state[i + 1]) << "\t\t"
           << elec_state[i] << '\n';
        it += 1;
    }
    to << line('-') << '\n';
}

void chem::print_geometry(std::ostream& to,
                          const std::vector<Element>& atoms,
                          const arma::mat& xyz,
                          const std::string& unit)
{
    chem::Format<char> line;
    chem::Format<double> fix;
    line.width(58).fill('-');
    fix.fixed().width(10);

    if (!atoms.empty()) {
        to << line('-') << '\n'
           << "Center\tAtomic\t\t    Coordinates/" << unit << '\n'
           << "Number\tSymbol\t   X\t\t   Y\t\t   Z\n"
           << line('-') << '\n';
        for (std::size_t i = 0; i < atoms.size(); ++i) {
            to << i + 1 << '\t' << atoms[i].atomic_symbol << '\t';
            for (arma::uword j = 0; j < xyz.n_cols; ++j) {
                to << fix(xyz(i, j)) << '\t';
            }
            to << '\n';
        }
        to << line('-') << '\n';
    }
}

void chem::print_atomic_masses(std::ostream& to,
                               const std::vector<Element>& atoms)
{
    chem::Format<double> fix;
    fix.fixed().width(10);

    chem::Format<int> gen;
    gen.width(3);

    double totmass = 0.0;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        to << "Center " << gen(i + 1) << " has atomic number "
           << gen(atoms[i].atomic_number) << " and mass "
           << fix(atoms[i].atomic_mass) << '\n';
        totmass += atoms[i].atomic_mass;
    }
    to << "Molecular mass:\t" << totmass << " amu\n";
}
