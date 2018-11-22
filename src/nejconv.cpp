////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2009-2018 Stig Rune Sellevag. All rights reserved.
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

#include <numlib/matrix.h>
#include <numlib/math.h>
#include <stdutils/stdutils.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

//------------------------------------------------------------------------------

struct Init_error : std::runtime_error {
    Init_error(const std::string& s) : std::runtime_error(s) {}
};

struct Conv_error : std::runtime_error {
    Conv_error(const std::string& s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

// Global declarations:

Numlib::Grid e_grid; // energy grid
Numlib::Grid j_grid; // angular momentum grid

std::vector<double> freq;     // vibrational frequencies in cm**-1
Numlib::Vec<double> vibr_nos; // vibrational number of states

//------------------------------------------------------------------------------

// Forward declarations:

// Initialize input data.
void init(const std::string& filename);

// Count number of vibrational states using Beyer-Swinehart algorithm.
void vibr_count();

// Convolute number of states.
void nej_conv(const std::string& filename);

//------------------------------------------------------------------------------

//
// Program for convolution of N(E,J) of transitional degrees of freedom
// with vibrational number of states.
//
int main(int argc, char* argv[])
{
    auto args = Stdutils::arguments(argc, argv);
    if (args.size() != 3) {
        std::cerr << "Usage: " << args[0] << " inp_file nej_file\n";
        return 1;
    }

    try {
        init(args[1]);
        vibr_count();
        nej_conv(args[2]);
    }
    catch (std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

//------------------------------------------------------------------------------

void init(const std::string& filename)
{
    std::ifstream from;
    Stdutils::fopen(from, filename);

    e_grid.set(from, "EnergyGrid");
    j_grid.set(from, "AngMomGrid");

    auto pos = Stdutils::find_token(from, "Frequencies");
    if (pos != -1) {
        double v;
        while (from >> v) {
            if (v > 0.0) {
                freq.push_back(v);
            }
            else {
                throw Init_error("bad frequency: " + std::to_string(v));
            }
        }
    }
    else {
        throw Init_error("could not find Frequencies in " + filename);
    }
}

void vibr_count()
{
    int nsize = e_grid.size();

    vibr_nos.resize(nsize);
    vibr_nos = 0.0;
    vibr_nos(0) = 1.0;

    for (auto fi : freq) {
        auto rj = Numlib::round<Index>(fi / e_grid.step());
        for (Index e = 0; e < vibr_nos.size() - rj; ++e) {
            vibr_nos(rj + e) += vibr_nos(e);
        }
    }
}

void nej_conv(const std::string& filename)
{
    std::ifstream from;
    Stdutils::fopen(from, filename);

    double ee = 0.0;
    double jj = 0.0;
    Numlib::Vec<double> nej(e_grid.size());

    for (Index j = 0; j < j_grid.size(); ++j) {
        for (Index e = 0; e < e_grid.size(); ++e) {
            from >> ee >> jj >> nej(e);
            if (!from) {
                throw Conv_error("cannot read N(E,J) data from " + filename);
            }
            if (ee != e_grid[e]) {
                throw Conv_error("bad E grid: " + std::to_string(ee) + ", " +
                                 std::to_string(e_grid[e]));
            }
            if (jj != j_grid[j]) {
                throw Conv_error("bad J grid: " + std::to_string(jj) + ", " +
                                 std::to_string(j_grid[j]));
            }
        }
        nej = Numlib::conv(vibr_nos, nej);
        for (Index e = 0; e < e_grid.size(); ++e) {
            std::cout << e_grid[e] << ' ' << jj << ' ' << nej(e) << '\n';
        }
    }
}

