// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/io.h>
#include <chem/periodic_table.h>
#include <numlib/constants.h>
#include <stdutils/stdutils.h>
#include <string>
#include <sstream>
#include <stdexcept>

void Chem::read_xyz_format(std::istream& from,
                           std::vector<Element>& atoms,
                           Numlib::Mat<double>& xyz,
                           std::string& title)
{
    // Get number of atoms:
    Index natoms;
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
    for (Index i = 0; i < natoms; ++i) {
        from >> symbol >> x >> y >> z;
        atoms[i] = Periodic_table::get_element(symbol);
        xyz(i, 0) = x;
        xyz(i, 1) = y;
        xyz(i, 2) = z;
    }
}

void Chem::read_zmat_format(std::istream& from,
                            std::vector<Element>& atoms,
                            Numlib::Vec<double>& distances,
                            Numlib::Vec<double>& angles,
                            Numlib::Vec<double>& dihedrals,
                            Numlib::Vec<Index>& bond_connect,
                            Numlib::Vec<Index>& angle_connect,
                            Numlib::Vec<Index>& dihedral_connect)
{
    atoms.clear();

    from.clear();
    from.seekg(0, std::ios_base::beg);

    std::string line;
    std::string symbol;

    Index natoms = 0;
    std::streamoff pos = 0;

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
        else if (Periodic_table::atomic_symbol_is_valid(symbol)) {
            atoms.push_back(Periodic_table::get_element(symbol));
            natoms += 1;
        }
        else if (symbol.empty()) {
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
    Index iat1;
    double distance;
    if (natoms > 1) {
        std::getline(from, line);
        std::istringstream iss(line);
        iss >> symbol >> iat1 >> distance;
        distances(1) = distance;
        bond_connect(1) = iat1 - 1;
    }
    Index iat2;
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
    Index iat3;
    double dihedral;
    if (natoms > 3) {
        for (Index i = 3; i < natoms; ++i) {
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

void Chem::read_mol_formula(std::istream& from,
                            std::vector<Mol_formula>& formula)
{
    std::string buf;
    from >> buf;

    auto n = std::stoi(buf);
    assert(n >= 1);
    formula.resize(n);

    char ch;
    std::string atom;
    int stoich;

    from >> ch;
    if (ch != '[') {
        throw std::runtime_error("'[' missing in molecular formula");
    }
    for (int i = 0; i < n; ++i) {
        from >> atom >> stoich >> ch;
        if (!from) {
            throw std::runtime_error("found no data for molecular formula");
        }
        if ((ch != ',') && (ch != ';') && (ch != ']')) {
            throw std::runtime_error("bad separator in molecular formula: " +
                                     std::to_string(ch));
        }
        if (ch == ']') {
            from.unget();
        }
        if (!Periodic_table::atomic_symbol_is_valid(atom)) {
            throw std::runtime_error("bad atomic symbol: " + atom);
        }
        if (stoich < 1) {
            throw std::runtime_error("bad stoichiometry: " +
                                     std::to_string(stoich));
        }
        formula[i].atom = atom;
        formula[i].stoich = stoich;
    }
    from >> ch;
    if (ch != ']') {
        throw std::runtime_error("']' missing in molecular formula");
    }
}

void Chem::print_xyz_format(std::ostream& to,
                            const std::vector<Element>& atoms,
                            const Numlib::Mat<double>& xyz,
                            const std::string& title)
{
    Stdutils::Format<double> fix;
    fix.fixed().width(10);

    to << atoms.size() << '\n';
    to << title << '\n';
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        to << atoms[i].atomic_symbol << '\t';
        for (Index j = 0; j < xyz.cols(); ++j) {
            to << fix(xyz(i, j)) << '\t';
        }
        to << '\n';
    }
}

void Chem::print_zmat_format(std::ostream& to,
                             const std::vector<Element>& atoms,
                             const Numlib::Vec<double>& distances,
                             const Numlib::Vec<double>& angles,
                             const Numlib::Vec<double>& dihedrals,
                             const Numlib::Vec<Index>& bond_connect,
                             const Numlib::Vec<Index>& angle_connect,
                             const Numlib::Vec<Index>& dihedral_connect)
{
    Stdutils::Format<double> dfix;
    dfix.fixed().width(10).precision(4);

    Stdutils::Format<Index> ifix;
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

void Chem::print_spin_orbit_states(std::ostream& to,
                                   const Numlib::Vec<int>& so_degen,
                                   const Numlib::Vec<double>& so_energy)
{
    Assert::dynamic(same_extents(so_degen, so_energy));

    Stdutils::Format<char> line;
    line.width(34).fill('-');

    Stdutils::Format<double> fix;
    fix.fixed().width(6).precision(2);

    to << "Spin-orbit states:\n"
       << line('-') << '\n'
       << " #\tEnergy/cm^-1\tDegeneracy\n"
       << line('-') << '\n';
    for (Index i = 0; i < so_degen.size(); ++i) {
        to << " " << i + 1 << '\t' << fix(so_energy(i)) << "\t\t" << so_degen(i)
           << '\n';
    }
    to << line('-') << '\n';
}

void Chem::print_geometry(std::ostream& to,
                          const std::vector<Element>& atoms,
                          const Numlib::Mat<double>& xyz,
                          const std::string& unit)
{
    Stdutils::Format<char> line;
    Stdutils::Format<double> fix;
    line.width(58).fill('-');
    fix.fixed().width(10);

    if (!atoms.empty()) {
        to << line('-') << '\n'
           << "Center\tAtomic\t\t    Coordinates/" << unit << '\n'
           << "Number\tSymbol\t   X\t\t   Y\t\t   Z\n"
           << line('-') << '\n';
        for (std::size_t i = 0; i < atoms.size(); ++i) {
            to << i + 1 << '\t' << atoms[i].atomic_symbol << '\t';
            for (Index j = 0; j < xyz.cols(); ++j) {
                to << fix(xyz(i, j)) << '\t';
            }
            to << '\n';
        }
        to << line('-') << '\n';
    }
}

void Chem::print_atomic_masses(std::ostream& to,
                               const std::vector<Element>& atoms)
{
    Stdutils::Format<double> fix;
    fix.fixed().width(10);

    Stdutils::Format<std::size_t> gen;
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

void Chem::print_center_of_mass(std::ostream& to,
                                const Numlib::Vec<double>& com)
{
    Stdutils::Format<double> fix;
    fix.fixed().width(8).precision(4);

    to << "Center of mass (X, Y, Z):  " << fix(com(0)) << ", " << fix(com(1))
       << ", " << fix(com(2)) << '\n';
}

void Chem::print_principal_moments(std::ostream& to,
                                   const Numlib::Vec<double>& pmom,
                                   const Numlib::Mat<double>& paxis)
{
    Stdutils::Format<char> line;
    line.width(54).fill('-');

    to << "\nPrincipal axes and moments of inertia in atomic units:\n"
       << line('-') << '\n';

    Stdutils::Format<double> fix;
    fix.fixed().width(12);

    to << "\t\tA\t     B\t\t  C\n"
       << "Eigenvalue: " << fix(pmom(0)) << ' ' << fix(pmom(1)) << ' '
       << fix(pmom(2)) << '\n';
    to << "     X      ";
    for (Index i = 0; i < pmom.size(); ++i) {
        to << fix(paxis(0, i)) << ' ';
    }
    to << "\n     Y      ";
    for (Index i = 0; i < pmom.size(); ++i) {
        to << fix(paxis(1, i)) << ' ';
    }
    to << "\n     Z      ";
    for (Index i = 0; i < pmom.size(); ++i) {
        to << fix(paxis(2, i)) << ' ';
    }
    to << '\n';
}

void Chem::print_rot_constants(std::ostream& to,
                               int sigma,
                               const std::string& symm,
                               const Numlib::Vec<double>& rotc)
{
    using namespace Numlib::Constants;

    const double gHz2icm = giga / (c_0 * 100.0);

    Stdutils::Format<char> line;
    line.width(21).fill('-');

    Stdutils::Format<double> fix;
    fix.fixed().width(14);

    if (rotc(0) > 0.0) {
        to << "\nRotational constants:\n" << line('-') << '\n';

        if (symm.find("linear") == 0) {
            to << fix(rotc(0)) << " GHz\n"
               << fix(rotc(0) * gHz2icm) << " cm^-1\n\n";
        }
        else {
            double ra = rotc(0);
            double rb = rotc(1);
            double rc = rotc(2);
            to << "\tA\t\tB\t\tC\n"
               << fix(ra) << '\t' << fix(rb) << '\t' << fix(rc) << " GHz\n"
               << fix(ra * gHz2icm) << '\t' << fix(rb * gHz2icm) << '\t'
               << fix(rc * gHz2icm) << " cm^-1\n\n";
        }
        to << "Rotational symmetry number: " << sigma << '\n';
        if (symm.find("atom") == 0) {
            to << "Rotational symmetry: This is an atom\n";
        }
        else {
            to << "Rotational symmetry: " << symm << '\n';
        }
    }
}
