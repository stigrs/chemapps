/**
   @file molrot.cpp

   This file is part of ChemApps - A C++ Chemistry Toolkit

   Copyright (C) 2009-2017  Stig Rune Sellevag

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

#include <chem/molecule.h>
#include <chem/utils.h>
#include <boost/program_options.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

/// Program for rotational analysis of molecules.
int main(int argc, char* argv[])
{
    namespace po = boost::program_options;

    po::options_description options("Allowed options");
    // clang-format off
    options.add_options()
        ("help,h", "display help message")
        ("file,f", po::value<std::string>(), "input file");
    // clang-format on
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    std::string input_file;

    if (vm.count("help")) {
        std::cout << options << '\n';
        return 0;
    }
    if (vm.count("file")) {
        input_file = vm["file"].as<std::string>();
    }
    else {
        std::cerr << options << '\n';
        return 1;
    }

    try {
        std::ifstream from;
        std::ofstream to;

        std::string output_file;
        output_file = chem::strip_suffix(input_file, ".inp");
        output_file = output_file + ".out";

        chem::fopen(from, input_file);
        chem::fopen(to, output_file.c_str());

        Molecule mol(from, to, "Molecule");
        mol.get_rot()->analysis(to);
        mol.get_tor()->analysis(to);
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}
