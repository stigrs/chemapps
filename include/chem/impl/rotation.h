// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

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

        Rotation(Geometry& g) : geom(g), sigma{1}, aligned{false} {}

        Rotation(std::istream& from, const std::string& key, Geometry& g);

        // Copy semantics:
        Rotation(const Rotation&) = default;
        Rotation& operator=(const Rotation&) = default;

        // Move semantics:
        Rotation(Rotation&&) = default;
        Rotation& operator=(Rotation&&) = default;

        ~Rotation() {}

        // Get rotational symmetry number.
        auto rot_sigma() const { return sigma; }

        // Compute rotational constants.
        Numlib::Vec<double> constants();

        // Compute rotational symmetry.
        std::string symmetry();

        // Get principal moments.
        auto principal_moments();

        // Get principal axes;
        auto principal_axes();

        // Rotate to principal axes.
        //
        // Note: This coordinate system is not the same as Gaussian's
        // standard orientation.
        //
        void rotate_to_principal_axes();

    private:
        // Move geometry to center of mass.
        void move_to_com();

        // Compute center of mass coordinates.
        Numlib::Vec<double> center_of_mass() const;

        // Compute principal moments of inertia.
        void calc_principal_moments();

        Geometry& geom;

        int sigma;
        bool aligned;

        Numlib::Vec<double> pmom = Numlib::zeros<Numlib::Vec<double>>(3);
        Numlib::Mat<double> paxis = Numlib::zeros<Numlib::Mat<double>>(3, 3);
    };

    inline void Rotation::move_to_com()
    {
        auto com = center_of_mass();
        Numlib::translate(geom.cart_coord(), -com(0), -com(1), -com(2));
    }

    inline void Rotation::rotate_to_principal_axes()
    {
        if (!geom.atoms().empty() && !aligned) {
            move_to_com();
            calc_principal_moments();
            aligned = true;
        }
    }

    inline auto Rotation::principal_moments()
    {
        if (!aligned) {
            rotate_to_principal_axes();
        }
        return pmom;
    }

    inline auto Rotation::principal_axes()
    {
        if (!aligned) {
            rotate_to_principal_axes();
        }
        return paxis;
    }

}  // namespace Impl

}  // namespace Chem

#endif  // CHEM_ROTATION_H
