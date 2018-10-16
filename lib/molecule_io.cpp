// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/impl/molecule_io.h>
#include <chem/periodic_table.h>
#include <stdutils/stdutils.h>

void Chem::Impl::read_xyz_format(std::istream& from,
                                 std::vector<Chem::Element>& atoms,
                                 Numlib::Mat<double>& xyz,
                                 std::string& title)
{
    // Get number of atoms:
    int natoms;
    from >> natoms;
    from.ignore(); // need to consume '\n' before reading title line
    from.clear();

    xyz.resize(natoms, 3);
    atoms.resize(natoms);

    std::string symbol;
    double x;
    double y;
    double z;

    // Read title line:
    std::getline(from, title);
    title = Stdutils::trim(title, " ");

    // Read XYZ coordinates:
    for (int i = 0; i < natoms; ++i) {
        from >> symbol >> x >> y >> z;
        atoms[i]  = Chem::Periodic_table::get_element(symbol);
        xyz(i, 0) = x;
        xyz(i, 1) = y;
        xyz(i, 2) = z;
    }
}

void Chem::Impl::read_zmat_format(std::istream& from,
                                  std::vector<Chem::Element>& atoms,
                                  Numlib::Vec<double>& distances,
                                  Numlib::Vec<double>& angles,
                                  Numlib::Vec<double>& dihedrals,
                                  Numlib::Vec<int>& bond_connect,
                                  Numlib::Vec<int>& angle_connect,
                                  Numlib::Vec<int>& dihedral_connect)
{
    atoms.clear();

    std::string line;
    std::string symbol;

    int natoms = 0;
    std::streamoff pos = 0;

    // Count number of atoms:

    while (std::getline(from, line)) {
        std::istringstream iss(line);
        iss >> symbol;
        if (symbol.find("#") != std::string::npos) {
            continue;
        }
        if (symbol.find("zmatrix") != std::string::npos) {
            pos = from.tellg();
        }
        if (Chem::Periodic_table::atomic_symbol_is_valid(symbol)) {
            atoms.push_back(Chem::Periodic_table::get_element(symbol));
            natoms += 1;
        }
        if (symbol.empty()) {
            break;
        }
    }
    // Load Z matrix:

    from.clear();
    from.seekg(pos, std::ios_base::beg);

    distances.resize(natoms);
    angles.resize(natoms);
    dihedrals.resize(natoms);
    bond_connect.resize(natoms);
    angle_connect.resize(natoms);
    dihedral_connect.resize(natoms);

    if (natoms > 0) {
        std::getline(from, line); // first atom is already read
    }
    int iat1;
    double distance;
    if (natoms > 1) {
        std::getline(from, line);
        std::istringstream iss(line);
        iss >> symbol >> iat1 >> distance;
        distances(1) = distance;
        bond_connect(1) = iat1 - 1;
    }
    int iat2;
    double angle;
    if (natoms > 2) {
        std::getline(from, line);
        std::istringstream iss(line);
        iss >> symbol >> iat1 >> distance >> iat2 >> angle;
        distances(2) = distance;
        bond_connect(2) = iat1 - 1;
        angles(2) = angle;
        angle_connect(2) = iat2 - 1;
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
            distances(i) = distance;
            bond_connect(i) = iat1 - 1;
            angles(i) = angle;
            angle_connect(i) = iat2 - 1;
            dihedrals(i) = dihedral;
            dihedral_connect(i) = iat3 - 1;
        }
    }
}

void Chem::Impl::print_elec_states(const Numlib::Vec<double>& elec_state)
{
    using namespace Stdutils;

    Format<char> line;
    line.width(34).fill('-');

    Format<double> fix;
    fix.fixed().width(6).precision(2);

    std::cout << "Electronic states:\n"
              << line('-') << '\n'
              << " #\tEnergy/cm^-1\tDegeneracy\n"
              << line('-') << '\n';
    int it = 1;
    for (Index i = 0; i < elec_state.size(); i += 2) {
        std::cout << " " << it << '\t' << fix(elec_state(i + 1)) << "\t\t"
                  << elec_state(i) << '\n';
        it += 1;
    }
    std::cout << line('-') << '\n';
}

void Chem::Impl::print_zmat_format(const std::vector<Chem::Element>& atoms,
                                   const Numlib::Vec<double>& distances,
                                   const Numlib::Vec<double>& angles,
                                   const Numlib::Vec<double>& dihedrals,
                                   const Numlib::Vec<int>& bond_connect,
                                   const Numlib::Vec<int>& angle_connect,
                                   const Numlib::Vec<int>& dihedral_connect)
{
    Stdutils::Format<double> dfix;
    dfix.fixed().width(10).precision(4);

    Stdutils::Format<int> ifix;
    ifix.fixed().width(3);

    if (!atoms.empty()) {
        std::cout << atoms[0].atomic_symbol << '\n';
    }
    if (atoms.size() > 1) {
        std::cout << atoms[1].atomic_symbol << '\t' << ifix(bond_connect(1) + 1)
                  << "  " << dfix(distances(1)) << '\n';
    }
    if (atoms.size() > 2) {
        std::cout << atoms[2].atomic_symbol << '\t' << ifix(bond_connect(2) + 1)
                  << "  " << dfix(distances(2)) << "  "
                  << ifix(angle_connect(2) + 1) << "  " << dfix(angles(2))
                  << '\n';
    }
    if (atoms.size() > 3) {
        for (std::size_t i = 3; i < atoms.size(); ++i) {
            std::cout << atoms[i].atomic_symbol << '\t'
                      << ifix(bond_connect(i) + 1) << "  " << dfix(distances(i))
                      << "  " << ifix(angle_connect(i) + 1) << "  "
                      << dfix(angles(i)) << "  "
                      << ifix(dihedral_connect(i) + 1) << "  "
                      << dfix(dihedrals(i)) << '\n';
        }
    }
}

