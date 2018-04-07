////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2018 Stig Rune Sellevag. All rights reserved.
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

#include <chem/gaussian.h>
#include <sstream>

Gauss_version Gaussian::get_version() const
{
    const char pattern_start[] = "Cite this work as:";
    const char pattern_ver[]   = "Gaussian";

    std::streampos orig_pos = from.tellg();

    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file
    from.clear();

    Gauss_version version = unknown;

    std::string line;
    std::string word;
    std::string ver;

    while (std::getline(from, line)) {
        if (line.find(pattern_start, 0) != std::string::npos) {
            std::getline(from, line);
            std::istringstream iss(line);
            iss >> word >> ver;
            if (word == pattern_ver) {
                if (ver == "94,") {
                    version = g94;
                }
                else if (ver == "98,") {
                    version = g98;
                }
                else if (ver == "03,") {
                    version = g03;
                }
                else if (ver == "09,") {
                    version = g09;
                }
            }
        }
    }
    if (version == unknown) {
        throw Gauss_error("unknown Gaussian version");
    }
    from.clear();
    from.seekg(orig_pos, std::ios_base::beg);  // move to original position
    from.clear();

    return version;
}

int Gaussian::get_natoms() const
{
    std::streampos orig_pos = from.tellg();

    from.clear();
    from.seekg(0, std::ios_base::beg);  // move to beginning of file
    from.clear();

    int natoms = 0;

    if (filetype == out) {
        const char pattern_start[] = "Z-Matrix";
        const char pattern_end[]   = "Distance";

        std::string line;
        std::string token;

        while (std::getline(from, line)) {
            std::istringstream iss1(line);
            iss1 >> token;
            if (token == pattern_start) {
                for (int i = 0; i < 4; ++i) {  // ignore four lines
                    from.ignore(256, '\n');
                }
                while (std::getline(from, line)) {
                    std::istringstream iss2(line);
                    if (line[1] == '-') {
                        break;
                    }
                    else {
                        iss2 >> natoms;
                    }
                }
            }
            else if (token == pattern_end) {
                break;
            }
        }
    }
    else {  // filetype == fchk
        const char pattern[] = "Number of atoms";

        std::string line;
        while (std::getline(from, line)) {
            std::string::size_type pos = line.find(pattern);
            if (pos != std::string::npos) {
                std::istringstream iss(line);
                char ignore;
                iss.ignore(std::strlen(pattern), '\n');
                iss >> ignore >> natoms;
                break;
            }
        }
    }
    if (natoms == 0) {
        throw Gauss_error(
            "could not determine number of atoms from Gaussian file");
    }
    from.clear();
    from.seekg(orig_pos, std::ios_base::beg);  // move to original position
    from.clear();

    return natoms;
}