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

#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <chem/molecule.h>
#include <chem/utils.h>
#include <chem/argparse.h>

/**
   Program for rotational analysis of molecules.
*/
int main(int argc, char* argv[]) 
{
    std::string input_file;

    Argparse::init(argc, argv);
    if (Argparse::has_switch("-f")) {
        input_file = Argparse::read("-f", "");
    } else if (Argparse::has_switch("-h")) {
        std::cout << "usage: " << argv[0] << " -f input_file.inp\n";
        return 0;
    } else {
        std::cerr << "usage: " << argv[0] << " -f input_file.inp\n";
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
