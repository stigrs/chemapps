// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_GEOMETRY_H
#define CHEM_GEOMETRY_H

#include <chem/element.h>
#include <chem/zmatrix.h>
#include <numlib/matrix.h>
#include <iostream>
#include <string>
#include <vector>

namespace Chem {

namespace Impl {

    // Class for handling molecular geometries.
    //
    class Geometry {
    public:
        Geometry() = default;

        Geometry(std::istream& from, const std::string& key);

        // Copy semantics:
        Geometry(const Geometry&) = default;
        Geometry& operator=(const Geometry&) = default;

        // Move semantics:
        Geometry(Geometry&&) = default;
        Geometry& operator=(Geometry&&) = default;

        ~Geometry() = default;

        // Get properties:

        double tot_mass() const;

        const auto& atoms() const { return atms; }

        auto& cart_coord() { return xyz; }
        const auto& cart_coord() const { return xyz; }

        auto& int_coord() { return zmat; }
        const auto& int_coord() const { return zmat; }

    private:
        std::vector<Chem::Element> atms;
        Numlib::Mat<double> xyz;
        Chem::Zmatrix zmat;
        std::string title;
    };

}  // namespace Impl

}  // namespace Chem

#endif  // CHEM_GEOMETRY_H
