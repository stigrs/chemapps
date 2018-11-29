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

#ifndef CHEM_ROTATION_H
#define CHEM_ROTATION_H

#include <chem/impl/geometry.h>
#include <numlib/math.h>
#include <numlib/matrix.h>
#include <iostream>
#include <string>

namespace Chem {

namespace Impl {

    // Class for handling molecular rotations.
    //
    class Rotation {
    public:
        Rotation() = delete;

        Rotation(Geometry& g)
            : geom(g), pmom(3), paxis(3, 3), aligned{false}, sigma{1}
        {
        }

        Rotation(std::istream& from, const std::string& key, Geometry& g);

        // Copy semantics:
        Rotation(const Rotation&) = default;
        Rotation& operator=(const Rotation&) = default;

        // Move semantics:
        Rotation(Rotation&&) = default;
        Rotation& operator=(Rotation&&) = default;

        ~Rotation() {}

        // Perform rotational analysis.
        void analysis(std::ostream& to) const;

        // Get rotational symmetry number.
        auto rot_sigma() const { return sigma; }

        // Compute rotational constants.
        Numlib::Vec<double> constants() const;

        // Compute rotational symmetry.
        std::string symmetry() const;

        // Get principal moments.
        auto principal_moments() const { return pmom; }

        // Get principal axes;
        auto principal_axes() const { return paxis; }

    private:
        // Move geometry to center of mass.
        void move_to_com();

        // Rotate to principal axes.
        //
        // Note: This coordinate system is not the same as Gaussian's
        // standard orientation.
        //
        void rotate_to_principal_axes();

        // Compute center of mass coordinates.
        Numlib::Vec<double> center_of_mass() const;

        // Compute principal moments of inertia.
        void calc_principal_moments();

        Geometry& geom;

        Numlib::Vec<double> pmom;
        Numlib::Mat<double> paxis;

        bool aligned;
        int sigma;
    };

    inline void Rotation::move_to_com()
    {
        auto com = center_of_mass();
        auto xyz = geom.cart_coord();
        Numlib::translate(xyz, -com(0), -com(1), -com(2));
        geom.set_cart_coord(xyz);
    }

    inline void Rotation::rotate_to_principal_axes()
    {
        if (!geom.atoms().empty() && !aligned) {
            move_to_com();
            calc_principal_moments();
            aligned = true;
        }
    }

} // namespace Impl

} // namespace Chem

#endif // CHEM_ROTATION_H
