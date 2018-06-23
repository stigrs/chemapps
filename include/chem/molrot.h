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

#ifndef CHEM_MOLROT_H
#define CHEM_MOLROT_H

#include <chem/element.h>
#include <srs/array.h>
#include <srs/math.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// Error reporting:

struct Molrot_error : std::runtime_error {
    Molrot_error(std::string s) : std::runtime_error(s) {}
};

// Forward declarations to allow friend declarations:

class Molvib;
class Torsion;

//
// Class for handling molecular rotations.
//
class Molrot {
public:
    Molrot(std::vector<Element>& atoms_, srs::dmatrix& xyz_)
        : atoms(atoms_), xyz(xyz_)
    {
    }

    Molrot(std::istream& from,
           const std::string& key,
           std::vector<Element>& atoms_,
           srs::dmatrix& xyz_)
        : atoms(atoms_), xyz(xyz_)
    {
        init(from, key);
    }

    Molrot(const Molrot& rot);

    // Perform rotational analysis.
    void analysis(std::ostream& to = std::cout);

    // Calculate total molecular mass.
    double tot_mass() const;

    // Compute rotational constants.
    srs::dvector constants();

    // Return rotational symmetry number.
    double get_sigma() const { return sigma; }

    // Return rotational symmetry.
    std::string symmetry();

    friend class Molvib;
    friend class Torsion;

protected:
    // Initialize.
    void init(std::istream& from, const std::string& key);

    // Move geometry to center of mass.
    void move_to_com();

    // Compute center of mass coordinates.
    srs::dvector center_of_mass() const;

    // Compute principal moments of inertia.
    void principal_moments();

    // Rotate to principal axes.
    //
    // Note: This coordinate system is not the same as Gaussian's
    // standard orientation.
    void rotate_to_principal_axes();

    // Print center of mass.
    void print_center_of_mass(std::ostream& to = std::cout) const;

    // Print principal moments and axes.
    void print_principal_moments(std::ostream& to = std::cout) const;

    // Print rotational constants.
    void print_constants(std::ostream& to = std::cout);

    std::vector<Element>& atoms;
    srs::dmatrix& xyz;

    srs::dvector pmom;
    srs::dmatrix paxis;

    double sigma;
    bool aligned = false;
};

inline void Molrot::move_to_com()
{
    srs::dvector com = center_of_mass();
    srs::translate(xyz, -com(0), -com(1), -com(2));
}

#endif  // CHEM_MOLROT_H
