//////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2017 Stig Rune Sellevag. All rights reserved.
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
///////////////////////////////////////////////////////////////////////////////

#include <chem/arma_io.h>
#include <chem/datum.h>
#include <chem/molvib.h>
#include <chem/utils.h>

Molvib::Molvib(std::istream& from, const std::string& key)
{
    bool found = chem::find_section(from, key);
    if (found) {
        std::string token;
        while (from >> token) {
            if (token == "End") {
                break;
            }
            else if (token == "vibr") {
                chem::read_vector(from, freqs);
            }
        }
    }
    else {
        throw Molvib_error("cannot find " + key + " section");
    }
}

void Molvib::print(std::ostream& to)
{
    chem::Format<char> line;
    line.width(26).fill('-');

    chem::Format<double> fix;
    fix.fixed().width(8).precision(2);

    if (freqs.size() > 0) {
        int it = 0;
        to << "\nVibrational modes (cm^-1):\n" << line('-') << '\n';
        for (arma::uword i = 0; i < freqs.size(); ++i) {
            to << fix(freqs(i));
            if ((it == 8) && (freqs.size() > 9)) {
                to << '\n';
            }
            it += 1;
        }
        double zpe = zero_point_energy();
        to << "\n\nZero-point vibrational energy: " << zpe / datum::au_to_icm
           << " Hartree\n";
    }
}
