////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2018 Stig Rune Sellevag. All rights reserved.
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

#ifndef CHEM_GAUSS_DATA_H
#define CHEM_GAUSS_DATA_H

#include <srs/array.h>
#include <srs/packed.h>
#include <iostream>
#include <string>
#include <vector>

//------------------------------------------------------------------------------

enum Gauss_filetype { out, fchk };
enum Gauss_version { unknown, g94, g98, g03, g09 };

//------------------------------------------------------------------------------

// Struct for holding coordinates.
struct Gauss_coord {
    int natoms;
    srs::ivector atnum;
    srs::dmatrix xyz;
};

// Struct for holding condensed NMR spectrum.
struct Gauss_NMR {
    srs::ivector number;
    std::string atom;
    srs::dvector shield;
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
    srs::dvector get_scf_zpe_energy() const;

    // Get optimized Cartesian coordinates.
    void get_opt_cart_coord(struct Gauss_coord& coord) const;

    // Get vibrational frequencies.
    void get_freqs(srs::dvector& freqs) const;

    // Get hessians.
    void get_hessians(srs::packed_dmatrix& hess) const;

    // Get data from relaxed PES scan.
    void get_pes_scan_data(std::string& scan_coord,
                           srs::dvector& coord,
                           srs::dvector& energy) const;

    // Get NMR data.
    void get_nmr_data(std::vector<Gauss_NMR>& nmr,
                      const std::string& nmr_method,
                      const double degen_tol) const;

    // Get number of calculated IRC points including starting point.
    int get_no_irc_points() const;

    // Get VMEP and SMEP values for each optimized IRC point.
    void get_irc_data(srs::dvector& mep) const;

    // Get geometries for each optimized IRC point.
    void get_irc_geom(srs::dvector& geom) const;

    // Get gradients for each optimized IRC point.
    void get_irc_grad(srs::dvector& grad) const;

    // Get Hessians for each optimized IRC point (only for output files).
    void get_irc_hess(srs::dvector& hess) const;

    // Get ModRedundant coordinate.
    std::string get_modredundant_coord() const;

    // Print optimized geometry.
    void print_opt_geom(std::ostream& to = std::cout) const;

private:
    std::istream& from;
    Gauss_filetype filetype;
};

#endif  // CHEM_GAUSS_DATA_H
