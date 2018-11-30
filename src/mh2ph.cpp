// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4018 4267) // caused by cxxopts.hpp
#endif

#include <stdutils/stdutils.h>
#include <cxxopts.hpp>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#ifdef _MSC_VER
#pragma warning(pop)
#endif

// Convert MOLPRO force constants to Polyrate format.
//
int main(int argc, char* argv[])
{
    // clang-format off
    cxxopts::Options options(argv[0], "Convert MOLPRO force constants to Polyrate format");
    options.add_options()
        ("h,help", "display help message")
        ("f,file", "input file", cxxopts::value<std::string>());
    // clang-format on

    auto args = options.parse(argc, argv);

    std::string input_file;

    if (args.count("help")) {
        std::cout << options.help({"", "Group"}) << '\n';
        return 0;
    }
    if (args.count("file")) {
        input_file = args["file"].as<std::string>();
    }
    else {
        std::cerr << options.help({"", "Group"}) << '\n';
        return 1;
    }

    try {
        std::ifstream from;
        Stdutils::fopen(from, input_file);

        std::vector<std::string> names;
        std::vector<std::string> lines;

        std::string token, data;
        while (from >> token) {
            std::getline(from, data);
            bool is_new = true;
            for (std::size_t i = 0; i < names.size(); ++i) {
                if (token == names[i]) {
                    is_new = false;
                    lines[i] += data + '\n';
                    break;
                }
            }
            if (is_new) {
                names.push_back(token);
                lines.push_back(data + '\n');
            }
        }
        for (auto l : lines) {
            std::cout << l;
        }
    }
    catch (std::exception& e) {
        std::cerr << "what: " << e.what() << '\n';
        return 1;
    }
}

