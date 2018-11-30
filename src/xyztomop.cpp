// Copyright (c) 2009-2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <stdutils/stdutils.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

//------------------------------------------------------------------------------

// Error reporting:

struct IO_error : std::runtime_error {
    IO_error(const std::string& s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

// Forward declarations:

void generate_mop(const std::string& tml_file, const std::string& xyz_file);
void get_xyz(const std::string& xyz_file);

//------------------------------------------------------------------------------

// Program for generating MOPAC input file from a template input file and by
// extracting Cartesian coordinates from a XYZ file.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 3) {
        std::cerr << "Usage: " << args[0] << " mopac.tml file.xyz\n\n"
                  << "mopac.tml: Template file for MOPAC input file\n"
                  << "file.xyz:  File with XYZ coordinates\n";
        return 1;
    }

    try {
        generate_mop(args[1], args[2]);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

// Generate MOPAC input file.
void generate_mop(const std::string& tml_file, const std::string& xyz_file)
{
    std::ifstream from(tml_file.c_str());
    if (!from) {
        throw IO_error("cannot open " + tml_file);
    }

    const std::string geom_here = "GEOMETRY_HERE";

    std::string line;
    bool found = false;

    while (std::getline(from, line)) {
        if (line.find(geom_here, 0) != std::string::npos) {
            get_xyz(xyz_file);
            found = true;
        }
        else {
            std::cout << line << '\n';
        }
    }
    if (!found) {
        throw IO_error("could not find keyword " + geom_here);
    }
}

// Extract Cartesian geometry from XYZ file.
void get_xyz(const std::string& xyz_file)
{
    std::ifstream from(xyz_file.c_str());
    if (!from) {
        throw IO_error("cannot open " + xyz_file);
    }

    int natoms;
    from >> natoms;
    if (natoms < 1) {
        throw IO_error("bad number of atoms: " + std::to_string(natoms));
    }

    std::string atom;
    double x;
    double y;
    double z;

    Stdutils::Format<double> fix8;
    fix8.fixed().width(15).precision(8);

    int iter = 0;
    while (from >> atom >> x >> y >> z) {
        ++iter;
        std::cout << "   " << atom << "   " << fix8(x) << " 1 " << fix8(y)
                  << " 1 " << fix8(z) << " 1\n";
    }

    if (iter != natoms) {
        throw IO_error("error reading XYZ coordinates");
    }
}

