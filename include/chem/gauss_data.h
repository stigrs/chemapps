// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef CHEM_GAUSS_DATA_H
#define CHEM_GAUSS_DATA_H

#include <numlib/matrix.h>
#include <iostream>
#include <string>
#include <vector>

namespace Chem {

//------------------------------------------------------------------------------

enum Gauss_filetype { out, fchk };
enum Gauss_version { unknown, g94, g98, g03, g09 };

//------------------------------------------------------------------------------

// Struct for holding coordinates.
struct Gauss_coord {
    int natoms;
    std::vector<int> atnum;
    Numlib::Mat<double> xyz;
};

// Struct for holding condensed NMR spectrum.
struct Gauss_NMR {
    std::vector<int> number;
    std::string atom;
    std::vector<double> shield;
};

// Struct for holding NMR data.
struct Gauss_NMR_line {
    int number;
    std::string atom;
    double shield;

    // Sort by shieldings.
    bool operator<(const Gauss_NMR_line& val) const
    {
        return shield < val.shield;
    }
};

//------------------------------------------------------------------------------

//
// Class providing methods for extracting data from Gaussian output files and
// formatted checkpoint files.
//
class Gauss_data {
public:
    Gauss_data(std::istream& from_, Gauss_filetype type)
        : from(from_), filetype(type)
    {
    }

    // Get Gaussian version number (only for output files).
    Gauss_version get_version() const;

    // Check termination.
    bool check_termination() const;

    // Check geometry convergence.
    bool check_opt_conv() const;

    // Get number of atoms.
    int get_natoms() const;

    // Get electronic and zero-point energies.
    std::vector<double> get_scf_zpe_energy() const;

    // Get optimized Cartesian coordinates.
    void get_opt_cart_coord(struct Gauss_coord& coord) const;

    // Get vibrational frequencies.
    void get_freqs(std::vector<double>& freqs) const;

    // Get hessians.
    void get_hessians(Numlib::Symm_mat<double, Numlib::lo>& hess) const;

    // Get data from relaxed PES scan.
    void get_pes_scan_data(std::string& scan_coord,
                           std::vector<double>& coord,
                           std::vector<double>& energy) const;

    // Get NMR data.
    void get_nmr_data(std::vector<Gauss_NMR>& nmr,
                      const std::string& nmr_method,
                      const double degen_tol) const;

    // Get number of calculated IRC points including starting point.
    int get_no_irc_points() const;

    // Get VMEP and SMEP values for each optimized IRC point.
    void get_irc_data(std::vector<double>& mep) const;

    // Get geometries for each optimized IRC point.
    void get_irc_geom(std::vector<double>& geom) const;

    // Get gradients for each optimized IRC point.
    void get_irc_grad(std::vector<double>& grad) const;

    // Get Hessians for each optimized IRC point (only for output files).
    void get_irc_hess(std::vector<double>& hess) const;

    // Get ModRedundant coordinate.
    std::string get_modredundant_coord() const;

    // Print optimized geometry.
    void print_opt_geom(std::ostream& to = std::cout) const;

private:
    std::istream& from;
    Gauss_filetype filetype;
};

} // namespace Chem

#endif // CHEM_GAUSS_DATA_H

