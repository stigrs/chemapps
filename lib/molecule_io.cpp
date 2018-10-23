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

#include <chem/impl/molecule_io.h>
#include <chem/periodic_table.h>
#include <stdutils/stdutils.h>
#include <sstream>
#include <stdexcept>

void Chem::Impl::read_xyz_format(std::istream& from,
                                 std::vector<Element>& atoms,
                                 Numlib::Mat<double>& xyz,
                                 std::string& title)
{
    // Get number of atoms:
    int natoms;
    from >> natoms;
    from.ignore();  // need to consume '\n' before reading title line
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
        atoms[i] = Periodic_table::get_element(symbol);
        xyz(i, 0) = x;
        xyz(i, 1) = y;
        xyz(i, 2) = z;
    }
}

void Chem::Impl::read_zmat_format(std::istream& from,
                                  std::vector<Element>& atoms,
                                  Numlib::Vec<double>& distances,
                                  Numlib::Vec<double>& angles,
                                  Numlib::Vec<double>& dihedrals,
                                  Numlib::Vec<int>& bond_connect,
                                  Numlib::Vec<int>& angle_connect,
                                  Numlib::Vec<int>& dihedral_connect)
{
    atoms.clear();

    from.clear();
    from.seekg(0, std::ios_base::beg);

    std::string line;
    std::string symbol;

    int natoms         = 0;
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

    distances.resize(natoms, 0.0);
    angles.resize(natoms, 0.0);
    dihedrals.resize(natoms, 0.0);
    bond_connect.resize(natoms, 0);
    angle_connect.resize(natoms, 0);
    dihedral_connect.resize(natoms, 0);

    if (natoms > 0) {
        std::getline(from, line);  // first atom is already read
    }
    int iat1;
    double distance;
    if (natoms > 1) {
        std::getline(from, line);
        std::istringstream iss(line);
        iss >> symbol >> iat1 >> distance;
        distances(1)    = distance;
        bond_connect(1) = iat1 - 1;
    }
    int iat2;
    double angle;
    if (natoms > 2) {
        std::getline(from, line);
        std::istringstream iss(line);
        iss >> symbol >> iat1 >> distance >> iat2 >> angle;
        distances(2)     = distance;
        bond_connect(2)  = iat1 - 1;
        angles(2)        = angle;
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
            distances(i)        = distance;
            bond_connect(i)     = iat1 - 1;
            angles(i)           = angle;
            angle_connect(i)    = iat2 - 1;
            dihedrals(i)        = dihedral;
            dihedral_connect(i) = iat3 - 1;
        }
    }
}

void Chem::Impl::read_mol_formula(std::istream& from,
                                  std::vector<Mol_formula>& formula)
{
    std::string buf;
    from >> buf;

    auto n = Stdutils::from_string<int>(buf);
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
                                     Stdutils::to_string(ch));
        }
        if (ch == ']') {
            from.unget();
        }
        if (!Periodic_table::atomic_symbol_is_valid(atom)) {
            throw std::runtime_error("bad atomic symbol: " + atom);
        }
        if (stoich < 1) {
            throw std::runtime_error("bad stoichiometry: " +
                                     Stdutils::to_string(stoich));
        }
        formula[i].atom   = atom;
        formula[i].stoich = stoich;
    }
    from >> ch;
    if (ch != ']') {
        throw std::runtime_error("']' missing in molecular formula");
    }
}

void Chem::Impl::print_xyz_format(std::ostream& to,
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

void Chem::Impl::print_zmat_format(std::ostream& to,
                                   const std::vector<Element>& atoms,
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

void Chem::Impl::print_elec_states(std::ostream& to,
                                   const Numlib::Vec<double>& elec_state)
{
    Stdutils::Format<char> line;
    line.width(34).fill('-');

    Stdutils::Format<double> fix;
    fix.fixed().width(6).precision(2);

    to << "Electronic states:\n"
       << line('-') << '\n'
       << " #\tEnergy/cm^-1\tDegeneracy\n"
       << line('-') << '\n';
    int it = 1;
    for (Index i = 0; i < elec_state.size(); i += 2) {
        to << " " << it << '\t' << fix(elec_state(i + 1)) << "\t\t"
           << elec_state(i) << '\n';
        it += 1;
    }
    to << line('-') << '\n';
}

void Chem::Impl::print_geometry(std::ostream& to,
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

void Chem::Impl::print_atomic_masses(std::ostream& to,
                                     const std::vector<Element>& atoms)
{
    Stdutils::Format<double> fix;
    fix.fixed().width(10);

    Stdutils::Format<int> gen;
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

