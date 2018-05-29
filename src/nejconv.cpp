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

#include <srs/array.h>
#include <srs/math.h>
#include <srs/utils.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>


//------------------------------------------------------------------------------

struct Init_error : std::runtime_error {
    Init_error(std::string s) : std::runtime_error(s) {}
};

struct Conv_error : std::runtime_error {
    Conv_error(std::string s) : std::runtime_error(s) {}
};

//------------------------------------------------------------------------------

// Global declarations:

Grid e_grid;  // energy grid
Grid j_grid;  // angular momentum grid

std::vector<double> freq;        // vibrational frequencies in cm**-1
srs::Array<double, 1> vibr_nos;  // vibrational number of states

//------------------------------------------------------------------------------

// Forward declarations:

// Initialize input data.
void init(const std::string& filename);

// Count number of vibrational states using Beyer-Swinehart algorithm.
void vibr_count();

// Convolute number of states.
void nej_conv(const std::string& filename);

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " inp_file nej_file\n";
        return 1;
    }

    try {
        init(argv[1]);
        vibr_count();
        nej_conv(argv[2]);
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
    srs::fopen(from, filename);

    e_grid.set(from, "EnergyGrid");
    j_grid.set(from, "AngMomGrid");

    if (srs::find_section(from, "Frequencies")) {
        double v;
        while (from >> v) {
            if (v > 0.0) {
                freq.push_back(v);
            }
            else {
                throw Init_error("bad frequency: " + srs::to_string(v));
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
    vibr_nos    = 0.0;
    vibr_nos[0] = 1.0;

    std::size_t rj;
    for (std::size_t i = 0; i < freq.size(); ++i) {
        rj = srs::nint(freq[i] / e_grid.step());
        for (std::size_t e = 0; e < vibr_nos.size() - rj; ++e) {
            vibr_nos[rj + e] += vibr_nos[e];
        }
    }
}

void nej_conv(const std::string& filename)
{
    std::ifstream from;
    srs::fopen(from, filename);

    double ee = 0.0;
    double jj = 0.0;
    srs::Array<double, 1> nej(e_grid.size());

    for (srs::size_t j = 0; j < j_grid.size(); ++j) {
        for (srs::size_t e = 0; e < e_grid.size(); ++e) {
            from >> ee >> jj >> nej[e];
            if (!from) {
                throw Conv_error("cannot read N(E,J) data from " + filename);
            }
            if (ee != e_grid[e]) {
                throw Conv_error("bad E grid: " + srs::to_string(ee) + ", "
                                 + srs::to_string(e_grid[e]));
            }
            if (jj != j_grid[j]) {
                throw Conv_error("bad J grid: " + srs::to_string(jj) + ", "
                                 + srs::to_string(j_grid[j]));
            }
        }
        nej = srs::conv(vibr_nos, nej);
        for (srs::size_t e = 0; e < e_grid.size(); ++e) {
            std::cout << e_grid[e] << ' ' << jj << ' ' << nej[e] << '\n';
        }
    }
}
