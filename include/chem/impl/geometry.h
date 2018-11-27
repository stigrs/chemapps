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

        std::string info() const { return title; }

        double tot_mass() const;

        const auto& atoms() const { return atms; }

        auto& cart_coord() { return xyz; }
        const auto& cart_coord() const { return xyz; }

        auto& int_coord() { return zmat; }
        const auto& int_coord() const { return zmat; }

        // Set properties:

        void set_cart_coord(const Numlib::Mat<double>& x);

    private:
        std::vector<Chem::Element> atms;
        Numlib::Mat<double> xyz;
        Chem::Zmatrix zmat;
        std::string title;
    };

    inline void Geometry::set_cart_coord(const Numlib::Mat<double>& x)
    {
        xyz = x;
        zmat.set(atms, xyz);
    }

} // namespace Impl

} // namespace Chem

#endif // CHEM_GEOMETRY_H
