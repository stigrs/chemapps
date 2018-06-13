////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////

#include <chem/molvib.h>
#include <srs/datum.h>
#include <srs/utils.h>


Molvib::Molvib(std::istream& from, const std::string& key)
{
    bool found = srs::find_section(from, key);
    if (found) {
        std::string token;
        while (from >> token) {
            if (token == "End") {
                break;
            }
            else if (token == "vibr") {
                from >> freqs;
            }
			else if (token == "hess") {
				from >> hess;
			}
        }
    }
    else {
        throw Molvib_error("cannot find " + key + " section");
    }
    // TODO (stigrs@gmail.com) Validate input data
}

void Molvib::print(std::ostream& to)
{
    srs::Format<char> line;
    line.width(26).fill('-');

    srs::Format<double> fix;
    fix.fixed().width(8).precision(2);

    if (freqs.size() > 0) {
        int it = 0;
        to << "\nVibrational modes (cm^-1):\n" << line('-') << '\n';
        for (srs::size_t i = 0; i < freqs.size(); ++i) {
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
