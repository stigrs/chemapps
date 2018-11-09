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
#include <chem/impl/io_support.h>
#include <stdutils/stdutils.h>

Chem::Molecule::Molecule(std::istream& from,
                         const std::string& key,
                         bool verbose)
    : elec(from, key),
      geom(from, key),
      rot(from, key, geom),
      vib(from, key, geom, rot),
      tor(from, key, geom, rot)
{
    if (verbose) {
        Stdutils::Format<char> line;
        line.width(15 + key.size()).fill('=');

        Stdutils::Format<double> fix;
        fix.fixed().precision(6);

        std::cout << "Input data on " << key << ":\n" << line('=') << '\n';
        std::cout << "Electronic energy: " << fix(elec.elec_energy())
                  << " Hartree\n"
                  << "Charge: " << elec.net_charge() << '\n'
                  << "Spin multiplicity: " << elec.spin_mult() << '\n';
        Chem::Impl::print_spin_orbit_states(std::cout, elec.spin_orbit_degen(),
                                            elec.spin_orbit_energy());
        std::cout << "\nInput orientation:\n";
        Chem::Impl::print_geometry(std::cout, geom.atoms(), geom.cart_coord());
        Chem::Impl::print_atomic_masses(std::cout, geom.atoms());
        vib.print(std::cout);
    }
}

