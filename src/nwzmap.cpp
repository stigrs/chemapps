////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2011-2018 Stig Rune Sellevag. All rights reserved.
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

#include <srs/utils.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

//------------------------------------------------------------------------------

// Forward declarations:

void get_map(std::istream& from, std::vector<int>& map);
double stretch(std::istream& from, int i, int j);
double bend(std::istream& from, int i, int j, int k);
double torsion(std::istream& from, int i, int j, int k, int l);
bool found(const std::vector<int>& a, std::vector<int>& b);

//------------------------------------------------------------------------------

//
// Program for mapping NWChem Z-matrix to different numbering of atoms.
//
int main(int argc, char* argv[])
{
    if (argc < 4) {
        std::cerr << "usage: " << argv[0] << " map_file tml_file zmat_file\n";
        return 1;
    }

    std::ifstream map_file(argv[1]);
    if (!map_file) {
        std::cerr << "could not open mapping file\n";
        return 1;
    }

    std::ifstream tml_file(argv[2]);
    if (!tml_file) {
        std::cerr << "could not open template file\n";
        return 1;
    }

    std::ifstream zmat_file(argv[3]);
    if (!zmat_file) {
        std::cerr << "could not open z-matrix file\n";
        return 1;
    }

    // Get mapping:

    std::vector<int> map;
    get_map(map_file, map);

    // Parse template file:

    unsigned iatom = 0;

    int j;
    int k;
    int l;

    double r;
    double a;
    double d;

    std::string atom;
    std::string rvar;
    std::string avar;
    std::string dvar;
    std::string line;

    srs::Format<double> fix5;
    fix5.fixed().precision(5);

    // Parse first atom:

    if (std::getline(tml_file, line)) {
        iatom++;
    }
    else {
        std::cerr << "template file is empty\n";
        return 1;
    }

    // Parse second atom:

    std::getline(tml_file, line);
    std::istringstream iss(line);
    if (iss >> atom >> j >> rvar) {
        iatom++;
        r = stretch(zmat_file, map[iatom - 1], map[j - 1]);
        std::cout << "  " << rvar << "  " << fix5(r) << '\n';
    }
    else {
        std::cerr << "bad template format, got " << line << '\n';
        return 1;
    }

    // Parse third atom:

    if (std::getline(tml_file, line)) {
        iss.clear();
        iss.str(line);
        if (iss >> atom >> j >> rvar >> k >> avar) {
            iatom++;
            r = stretch(zmat_file, map[iatom - 1], map[j - 1]);
            a = bend(zmat_file, map[iatom - 1], map[j - 1], map[k - 1]);
            std::cout << "  " << rvar << "  " << fix5(r) << '\n'
                      << "  " << avar << "  " << fix5(a) << '\n';
        }
        else {
            std::cerr << "bad template format, got " << line << '\n';
            return 1;
        }
    }

    // Parse the rest of the atoms:

    while (std::getline(tml_file, line)) {
        iss.clear();
        iss.str(line);
        if (iss >> atom >> j >> rvar >> k >> avar >> l >> dvar) {
            iatom++;
            r = stretch(zmat_file, map[iatom - 1], map[j - 1]);
            a = bend(zmat_file, map[iatom - 1], map[j - 1], map[k - 1]);
            d = torsion(
                zmat_file, map[iatom - 1], map[j - 1], map[k - 1], map[l - 1]);
            std::cout << "  " << rvar << "  " << fix5(r) << '\n'
                      << "  " << avar << "  " << fix5(a) << '\n'
                      << "  " << dvar << "  " << fix5(d) << '\n';
        }
        else {
            std::cerr << "bad template format, got " << line << '\n';
            return 1;
        }
    }
    if (iatom != map.size()) {
        std::cerr << "incompatible number of atoms in map and template file\n";
        return 1;
    }
}

//------------------------------------------------------------------------------

void get_map(std::istream& from, std::vector<int>& map)
{
    int iter;
    while (from >> iter) {
        map.push_back(iter);
    }
}

double stretch(std::istream& from, int i, int j)
{
    from.clear();
    from.seekg(0, std::ios_base::beg);

    std::vector<int> va(2);
    va[0] = i;
    va[1] = j;
    std::sort(va.begin(), va.end());

    int it;
    int ii;
    int jj;

    double value;

    std::string token;
    std::string line;

    std::vector<int> vb(2);

    while (std::getline(from, line)) {
        std::istringstream iss(line);
        while (iss >> it >> token) {
            if (token == "Stretch") {
                while (iss >> ii >> jj >> value) {
                    vb[0] = ii;
                    vb[1] = jj;
                    std::sort(vb.begin(), vb.end());
                    if (found(va, vb)) {
                        return value;
                    }
                }
            }
        }
    }
    return 0.0;
}

double bend(std::istream& from, int i, int j, int k)
{
    from.clear();
    from.seekg(0, std::ios_base::beg);

    std::vector<int> va(3);
    va[0] = i;
    va[1] = j;
    va[2] = k;
    std::sort(va.begin(), va.end());

    int it;
    int ii;
    int jj;
    int kk;

    double value;

    std::string token;
    std::string line;

    std::vector<int> vb(3);

    while (std::getline(from, line)) {
        std::istringstream iss(line);
        while (iss >> it >> token) {
            if (token == "Bend") {
                while (iss >> ii >> jj >> kk >> value) {
                    vb[0] = ii;
                    vb[1] = jj;
                    vb[2] = kk;
                    std::sort(vb.begin(), vb.end());
                    if (found(va, vb)) {
                        return value;
                    }
                }
            }
        }
    }
    return 0.0;
}

double torsion(std::istream& from, int i, int j, int k, int l)
{
    from.clear();
    from.seekg(0, std::ios_base::beg);

    std::vector<int> va(4);
    va[0] = i;
    va[1] = j;
    va[2] = k;
    va[3] = l;
    std::sort(va.begin(), va.end());

    int it;
    int ii;
    int jj;
    int kk;
    int ll;

    double value;

    std::string token;
    std::string line;

    std::vector<int> vb(4);

    while (std::getline(from, line)) {
        std::istringstream iss(line);
        while (iss >> it >> token) {
            if (token == "Torsion") {
                while (iss >> ii >> jj >> kk >> ll >> value) {
                    vb[0] = ii;
                    vb[1] = jj;
                    vb[2] = kk;
                    vb[3] = ll;
                    std::sort(vb.begin(), vb.end());
                    if (found(va, vb)) {
                        return value;
                    }
                }
            }
        }
    }
    return 0.0;
}

bool found(const std::vector<int>& a, std::vector<int>& b)
{
    std::sort(b.begin(), b.end());
    return std::equal(a.begin(), a.end(), b.begin());
}
