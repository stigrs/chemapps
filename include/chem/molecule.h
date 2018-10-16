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

#ifndef CHEM_MOLECULE_H
#define CHEM_MOLECULE_H

#include <chem/impl/elec_state.h>
#include <chem/impl/geometry.h>
#include <chem/impl/rotation.h>
#include <chem/impl/vibration.h>
#include <iostream>
#include <string>

namespace Chem {

// Class for holding molecule objects.
//
class Molecule {
public:
    Molecule() = default;

    Molecule(std::istream& from,
             const std::string& key = "Molecule",
             bool verbose = false);

    // Copy semantics:
    Molecule(const Molecule&) = default;
    Molecule& operator=(const Molecule&) = default;

    // Move semantics:
    Molecule(Molecule&&) = default;
    Molecule& operator=(Molecule&&) = default;

    ~Molecule() = default;

    // Get number of atoms.
    auto num_atoms() const { return geom.atoms().size(); }

    // Get net electronic charge.
    auto net_charge() const { return elec.net_charge(); }

    // Get spin multiplicity.
    auto spin_mult() const { return elec.spin_mult(); }

    // Get electronic energy.
    auto elec_energy() const { return elec.elec_energy(); }

    // Get degeneracies of spin-orbit states.
    const auto& spin_orbit_degen() const { return elec.spin_orbit_degen(); }

    // Get energies of spin-orbit states.
    const auto& spin_orbit_energy() const { return elec.spin_orbit_energy(); }

    // Get atoms in molecule.
    const auto& atoms() const { return geom.atoms(); }

    // Get Cartesian coordinates.
    auto& cart_coord() { return geom.cart_coord(); }
    const auto& cart_coord() const { return geom.cart_coord(); }

    // Get internal coordinates.
    auto& int_coord() { return geom.int_coord(); }
    const auto& int_coord() const { return geom.int_coord(); }

    // Get rotational constants.
    auto rot_const() { return rot.constants(); }

    // Get rotational symmetry.
    auto rot_symmetry() { return rot.symmetry(); }

    // Get principal moments.
    auto principal_moments() { return rot.principal_moments(); }

    // Get principal axes.
    auto principal_axes() { return rot.principal_axes(); }

    // Get Hessians.
    const auto& hessians() const { return vib.hessians(); }

    // Zero-point vibrational energy.
    double zero_point_energy() const { return vib.zero_point_energy(); }

    // Vibrational frequencies.
    const auto& frequencies() const { return vib.frequencies(); }

    // Reduced masses of vibrational modes.
    const auto& vib_red_masses() const { return vib.red_masses(); }

    // Force constants of vibrational modes.
    const auto& vib_force_constants() const { return vib.force_constants(); }

private:
    Impl::Elec_state elec;  // molecular electronic states
    Impl::Geometry geom;    // molecular geometry
    Impl::Rotation rot;     // molecular rotations
    Impl::Vibration vib;    // molecular vibrations
};

}  // namespace Chem

#endif  // CHEM_MOLECULE_H

