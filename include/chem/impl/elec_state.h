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

#ifndef CHEM_ELEC_STATE_H
#define CHEM_ELEC_STATE_H

#include <numlib/matrix.h>
#include <iostream>
#include <string>

namespace Chem {

namespace Impl {

    // Class for handling electronic states.
    //
    class Elec_state {
    public:
        Elec_state()
            : charge{0}, spin{1}, energy{0.0}, so_degen{1}, so_energy{0.0}
        {
        }

        Elec_state(std::istream& from, const std::string& key);

        // Copy semantics:
        Elec_state(const Elec_state&) = default;
        Elec_state& operator=(const Elec_state&) = default;

        // Move semantics:
        Elec_state(Elec_state&&) = default;
        Elec_state& operator=(Elec_state&&) = default;

        ~Elec_state() = default;

        // Get properties:

        auto net_charge() const { return charge; }
        auto spin_mult() const { return spin; }
        auto elec_energy() const { return energy; }

        const auto& spin_orbit_degen() const { return so_degen; }
        const auto& spin_orbit_energy() const { return so_energy; }

        // Set properties:

        void set_net_charge(int value) { charge = value; }
        void set_spin_mult(int value) { spin = value; }
        void set_elec_energy(double value) { energy = value; }

    private:
        int charge;    // net electronic charge
        int spin;      // spin multiplicity
        double energy; // electronic energy

        Numlib::Vec<int> so_degen;     // degeneracies of spin-orbit states
        Numlib::Vec<double> so_energy; // energies of spin-orbit states
    };

} // namespace Impl

} // namespace Chem

#endif // CHEM_ELEC_STATE_H
