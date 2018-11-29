// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_VIBRATION_H
#define CHEM_VIBRATION_H

#include <chem/element.h>
#include <numlib/math.h>
#include <numlib/matrix.h>
#include <stdutils/stdutils.h>
#include <iostream>
#include <string>
#include <vector>

namespace Chem {

// Class for handling molecular vibrations.
//
class Vibration {
public:
    Vibration() = default;

    Vibration(const Numlib::Vec<double>& f) : freqs(f) {}

    Vibration(std::istream& from,
              const std::string& key,
              const std::vector<Element>& at,
              const Numlib::Mat<double>& x,
              const Numlib::Mat<double>& p);

    // Copy semantics:
    Vibration(const Vibration&) = default;
    Vibration& operator=(const Vibration&) = default;

    // Move semantics:
    Vibration(Vibration&&) = default;
    Vibration& operator=(Vibration&&) = default;

    ~Vibration() = default;

    // Perform vibrational analysis.
    void analysis(std::ostream& to = std::cout) const;

    // Get Hessians.
    const auto& hessians() const { return hess; }

    // Get vibrational frequencies.
    const auto& frequencies() const { return freqs; }

    // Calculate zero-point vibrational energy.
    double zero_point_energy() const;

    // Reduced masses for vibrational modes.
    const auto& red_masses() const { return mu_freqs; }

    // Force constants for vibrational modes.
    const auto& force_constants() const { return k_fc; }

    // Print vibrational modes.
    void print(std::ostream& to) const;

private:
    // Calculate mass-weighted Hessians.
    Numlib::Mat<double> mw_hessians() const;

    // Calculate vibrational normal modes.
    void calc_normal_modes();

    // Set up coodinate vectors for translation and rotation about
    // principal axes of inertia.
    void trans_rot_vec(Numlib::Cube<double>& dmat, int& n_tr_rot) const;

    // Transform Cartesian Hessians to internal coordinates.
    void trans_hess_int_coord(Numlib::Cube<double>& dmat,
                              Numlib::Mat<double>& lmat,
                              int& n_tr_rot) const;

    // Shuffle n_vib orthogonal vectors to the beginning of D matrix.
    void shuffle(Numlib::Cube<double>& dmat, int n_tr_rot) const;

    // Convert frequencies from atomic units to cm^-1.
    void freqs_unit_conv(Numlib::Vec<double>& vib) const;

    // Print Cartesian frequencies.
    void print_cart_freqs(std::ostream& to) const;

    // Print normal modes.
    void print_normal_modes(std::ostream& to) const;

    std::vector<Element> atms;
    Numlib::Mat<double> xyz;
    Numlib::Mat<double> paxis;

    Numlib::Symm_mat<double, Numlib::lo> hess; // packed Hessians
    Numlib::Vec<double> freqs;                 // vibrational frequencies
    Numlib::Vec<double> mu_freqs; // reduces masses for vibrational modes
    Numlib::Vec<double> k_fc;     // force constants for vibrational modes
    Numlib::Cube<double> l_cart;  // Cartesian displacements
};

inline void Vibration::analysis(std::ostream& to) const
{
    Stdutils::Format<char> line;
    line.width(21).fill('=');

    to << "\nVibrational analysis:\n" << line('=') << "\n\n";
    print(to);
}

} // namespace Chem

#endif // CHEM_VIBRATION_H

