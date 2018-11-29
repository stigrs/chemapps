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

#include <chem/molecule.h>
#include <chem/io.h>
#include <stdutils/stdutils.h>

Chem::Molecule::Molecule(std::istream& from,
                         std::ostream& to,
                         const std::string& key,
                         bool verbose)
    : geom_(from, key)
{
    if (verbose) {
        Stdutils::Format<char> line;
        line.width(15 + key.size()).fill('=');

        Stdutils::Format<double> fix;
        fix.fixed().precision(6);

        to << "Input data on " << key << ":\n" << line('=') << '\n';
#if 0
        to << "Electronic energy: " << fix(elec.elec_energy()) << " Hartree\n"
           << "Charge: " << elec.net_charge() << '\n'
           << "Spin multiplicity: " << elec.spin_mult() << '\n';
        Chem::Impl::print_spin_orbit_states(to, elec.spin_orbit_degen(),
                                            elec.spin_orbit_energy());
#endif
        to << "\nInput orientation:\n";
        Chem::print_geometry(to, geom_.atoms(), geom_.get_xyz());
        Chem::print_atomic_masses(to, geom_.atoms());
        // vib.print(to);
    }
}

