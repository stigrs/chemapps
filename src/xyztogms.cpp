// Copyright (c) 2009-2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/periodic_table.h>
#include <stdutils/stdutils.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

//------------------------------------------------------------------------------

// Forward declarations:

void convert(const std::string& input_file);

//------------------------------------------------------------------------------

// Converts a XYZ file to GAMESS $DATA format.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 2) {
        std::cerr << "usage: " << args[0] << " file.xyz\n";
        return 1;
    }

    try {
        convert(args[1]);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
    }
}

//------------------------------------------------------------------------------

void convert(const std::string& input_file)
{
    std::ifstream from;
    Stdutils::fopen(from, input_file);

    std::string line;
    std::string atsymb;
    double x;
    double y;
    double z;

    namespace Pt = Chem::Periodic_table;

    while (from >> atsymb >> x >> y >> z) {
        Stdutils::Format<double> fix(1);
        fix.fixed();
        double atnum = Pt::get_atomic_number(atsymb);
        std::cout << atsymb << "  " << fix(atnum) << "  " << x << "  " << y
                  << "  " << z << '\n';
    }
}

