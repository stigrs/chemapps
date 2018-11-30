// Copyright (c) 2009-2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <stdutils/stdutils.h>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

//------------------------------------------------------------------------------

struct Error : std::runtime_error {
    Error(const std::string& s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

bool ts_is_linear();
int get_ts_nmodes();
void punch_path_data();

std::ifstream from;

//------------------------------------------------------------------------------

// Prepares Polyrate path data for plotting.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 2) {
        std::cout << "usage: " << args[0] << " poly.fu6\n";
        return 1;
    }

    from.open(args[1]);
    if (!from.is_open()) {
        std::cerr << "cannot open " << args[1] << '\n';
        return 1;
    }

    try {
        punch_path_data();
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

bool ts_is_linear()
{
    from.clear();
    from.seekg(0, std::ios_base::beg); // move to beginning of stream

    bool linear = false;

    std::string line, word, tmp1, tmp2, tmp3;

    while (std::getline(from, line)) {
        if (line.find("Starting point Parameters:", 0) != std::string::npos) {
            while (std::getline(from, line)) {
                if (line.find("SPECIES:  species type", 0) !=
                    std::string::npos) {
                    std::istringstream iss(line);
                    iss >> tmp1 >> tmp2 >> tmp3 >> word;
                    if (word == "nonlints") {
                        linear = false;
                    }
                    else {
                        linear = true;
                    }
                    break;
                }
            }
        }
    }
    return linear;
}

int get_ts_nmodes()
{
    from.clear();
    from.seekg(0, std::ios_base::beg); // move to beginning of file

    int nmodes = 0;

    std::string line, word, ignore;

    while (std::getline(from, line)) {
        if (line.find("* Saddle point *", 0) != std::string::npos) {
            while (from >> word) {
                if (word == "NDIM") {
                    from >> ignore >> nmodes;
                    break;
                }
            }
        }
    }

    if (ts_is_linear()) {
        nmodes -= 6;
    }
    else {
        nmodes -= 7;
    }
    if (nmodes < 1) { // something is wrong...
        throw Error("bad number of generalized normal modes: " +
                    std::to_string(nmodes));
    }
    return nmodes;
}

void punch_path_data()
{
    int nmodes = get_ts_nmodes();

    from.clear();
    from.seekg(0, std::ios_base::beg); // move to beginning of stream

    std::string line, word;
    float smep, vmep, vag, dummy;
    std::vector<std::string> freq(nmodes); // use string to handle imag. freqs.

    while (std::getline(from, line)) {
        if (line.find("Classical and adiabatic energies", 0) !=
            std::string::npos) {
            std::getline(from, line); // ignore two lines
            std::getline(from, line);
            std::cout << "# " << line << '\n';
            while (from >> smep >> vmep >> vag >> dummy) {
                for (int i = 0; i < nmodes; i++) {
                    from >> freq[i];
                }
                std::cout << smep << '\t' << vmep << '\t' << vag << '\t'
                          << dummy;
                for (int i = 0; i < nmodes; i++) {
                    word = freq[i];
                    // write imaginary frequencies as negative frequencies
                    if (word.find('i', 0) != std::string::npos) {
                        std::cout << "\t-";
                        for (std::size_t j = 0; j < word.size() - 1; j++) {
                            std::cout << word[j];
                        }
                    }
                    else {
                        std::cout << '\t' << word;
                    }
                }
                std::cout << '\n';
            }
            from.clear(); // there could be several path data sections...
            std::cout << '\n';
        }
    }
}

