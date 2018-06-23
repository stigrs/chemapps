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

#include <chem/molvib.h>
#include <mkl.h>
#include <srs/datum.h>
#include <srs/math.h>
#include <srs/utils.h>
#include <cmath>


void Molvib::analysis(std::ostream& to)
{
    srs::Format<char> line;
    line.width(21).fill('=');

    to << "\nVibrational analysis:\n" << line('=') << "\n\n";

    if (hess.empty()) {
        print(to);
    }
    else {
        print_cart_freq(to);
        print_norm_modes(to);
    }
}

srs::packed_dmatrix Molvib::get_mw_hessians() const
{
    srs::packed_dmatrix hess_mw(hess);
    if (!hess_mw.empty()) {
        for (srs::size_t j = 0; j < hess_mw.cols(); ++j) {
            for (srs::size_t i = 0; i < j + 1; ++i) {
                hess_mw(i, j) /= std::sqrt(rot.atoms[i / 3].atomic_mass
                                           * rot.atoms[j / 3].atomic_mass);
            }
        }
    }
    return hess_mw;
}

srs::dvector Molvib::calc_cart_freqs() const
{
    srs::packed_dmatrix hess_mw = get_mw_hessians();

    srs::dvector freqs_cart(hess_mw.rows());
    srs::dmatrix v;

    srs::eigs(hess_mw, v, freqs_cart);
    freqs_unit_conv(freqs_cart);  // atomic units to cm^-1

    return freqs_cart;
}

double Molvib::zero_point_energy() const
{
    double zpe = 0.0;
    for (auto v : freqs) {
        if (v < 0.0) {  // ignore imaginary frequencies
            continue;
        }
        else {
            zpe += v;
        }
    }
    return 0.5 * zpe;
}

void Molvib::norm_mode_analysis()
{
    // Transform Cartesian force constants to internal coordinates:

    srs::dcube dmat;       // D matrix
    srs::dmatrix lmat;     // L matrix
    srs::size_t n_tr_rot;  // number of translational and rotational modes


    trans_hess_int_coord(dmat, lmat, n_tr_rot);
    srs::packed_dmatrix fc_int(lmat);

    // Calculate vibrational frequencies:

    srs::size_t natoms  = rot.atoms.size();
    srs::size_t natoms3 = 3 * natoms;
    srs::size_t n_vib   = natoms3 - n_tr_rot;

    freqs.resize(n_vib);
    lmat.resize(n_vib, n_vib);

    srs::eigs(fc_int, lmat, freqs);
    freqs_unit_conv(freqs);

    // Calculate reduced masses, force constants and Cartesian displacements:

    const double factor = 4.0 * std::pow(datum::c_0 * 100.0, 2.0) * datum::m_u
                          * datum::pi * datum::pi / 100.0;

    mu_freqs.resize(n_vib);
    k_fc.resize(n_vib);

    srs::dmatrix tmp(3, natoms);
    srs::dmatrix ltmp(natoms3, n_vib);

    for (srs::size_t i = 0; i < n_vib; ++i) {
        double y = 0.0;
        for (srs::size_t k1 = 0; k1 < natoms; ++k1) {
            for (srs::size_t k2 = 0; k2 < 3; ++k2) {
                double x = 0.0;
                for (srs::size_t j = 0; j < n_vib; ++j) {
                    x += dmat(k2, k1, j) * lmat(j, i);
                }
                x /= std::sqrt(rot.atoms[k1].atomic_mass);
                y += x * x;
                tmp(k2, k1) = x;
            }
        }
        mu_freqs(i) = 1.0 / y;
        k_fc(i)     = freqs(i) * freqs(i) * mu_freqs(i) * factor;
        y           = 1.0 / std::sqrt(y);
        for (srs::size_t it = 0; it < natoms3; ++it) {
            ltmp(it, i) = y * tmp.data()[it];
        }
    }
    l_cart = srs::dcube(3, natoms, n_vib, ltmp.data());
}

void Molvib::print(std::ostream& to)
{
    srs::Format<char> line;
    line.width(26).fill('-');

    srs::Format<double> fix;
    fix.fixed().width(8).precision(2);

    if (freqs.size() > 0) {
        int it = 0;
        to << "\nVibrational modes (cm^-1):\n" << line('-') << '\n';
        for (srs::size_t i = 0; i < freqs.size(); ++i) {
            to << fix(freqs(i));
            if ((it == 8) && (freqs.size() > 9)) {
                to << '\n';
            }
            it += 1;
        }
        double zpe = zero_point_energy();
        to << "\n\nZero-point vibrational energy: " << zpe / datum::au_to_icm
           << " Hartree\n";
    }
}

//------------------------------------------------------------------------------

void Molvib::init(std::istream& from, const std::string& key)
{
    bool found = srs::find_section(from, key);
    if (found) {
        std::string token;
        srs::dvector tmp;
        while (from >> token) {
            if (token == "End") {
                break;
            }
            else if (token == "vibr") {
                from >> freqs;
            }
            else if (token == "hessians") {
                from >> tmp;
                hess = srs::packed_dmatrix(tmp);
            }
        }
    }
    else {
        throw Molvib_error("cannot find " + key + " section");
    }
    // TODO (stigrs@gmail.com) Validate input data

    if (!hess.empty()) {
        norm_mode_analysis();
    }
}

void Molvib::trans_rot_vec(srs::dcube& dmat, int& n_tr_rot) const
{
    // This function generates the vectors corresponding to translations and
    // infinitesimal rotations. The variable n_tr_rot is set to 5 or 6,
    // depending on how many such vectors there are. The function does not
    // take the symmetry of the molecule into account.

    srs::size_t natoms = rot.atoms.size();
    srs::size_t n3     = 3 * natoms;

    // Determine center of mass, moments of inertia and rotation generators:

    rot.rotate_to_principal_axes();

    auto xyz   = rot.xyz * (1.0 / datum::a_0);
    auto paxis = rot.paxis;

    // Generate transformation matrix:

    for (srs::size_t i = 0; i < natoms; ++i) {
        double m      = std::sqrt(rot.atoms[i].atomic_mass);
        double cx     = xyz(i, 0);
        double cy     = xyz(i, 1);
        double cz     = xyz(i, 2);
        double cxp    = cx * paxis(0, 0) + cy * paxis(1, 0) + cz * paxis(2, 0);
        double cyp    = cx * paxis(0, 1) + cy * paxis(1, 1) + cz * paxis(2, 1);
        double czp    = cx * paxis(0, 2) + cy * paxis(1, 2) + cz * paxis(2, 2);
        dmat(0, i, 0) = m;
        dmat(1, i, 1) = m;
        dmat(2, i, 2) = m;
        dmat(0, i, 3) = (cyp * paxis(0, 2) - czp * paxis(0, 1)) * m;
        dmat(1, i, 3) = (cyp * paxis(1, 2) - czp * paxis(1, 1)) * m;
        dmat(2, i, 3) = (cyp * paxis(2, 2) - czp * paxis(2, 1)) * m;
        dmat(0, i, 4) = (czp * paxis(0, 0) - cxp * paxis(0, 2)) * m;
        dmat(1, i, 4) = (czp * paxis(1, 0) - cxp * paxis(1, 2)) * m;
        dmat(2, i, 4) = (czp * paxis(2, 0) - cxp * paxis(2, 2)) * m;
        dmat(0, i, 5) = (cxp * paxis(0, 1) - cyp * paxis(0, 0)) * m;
        dmat(1, i, 5) = (cxp * paxis(1, 1) - cyp * paxis(1, 0)) * m;
        dmat(2, i, 5) = (cxp * paxis(2, 1) - cyp * paxis(2, 0)) * m;
    }

    const double cutoff = 1.0e-12;

    n_tr_rot = 6;
    int i    = 1;
    while (i < n_tr_rot) {
        auto di  = dmat.depth(i - 1).flatten();
        double x = srs::dot(di, di);
        if (x < cutoff) {
            --n_tr_rot;
            auto tmp = dmat.depth(n_tr_rot).flatten();
            srs::copy_n(n3, tmp, di);
        }
        else {
            di *= (1.0 / std::sqrt(x));
            ++i;
        }
    }
    if (natoms == 1) {
        if (n_tr_rot != 3) {
            throw Molvib_error("single atom, but n_tr_rot is not 3");
        }
    }
}

void Molvib::trans_hess_int_coord(srs::dcube& dmat,
                                  srs::dmatrix& lmat,
                                  srs::size_t& n_tr_rot) const
{
    // Set up coordinate vectors for translation and rotation about principal
    // axes of inertia.

    // Since these coordinate vectors involve mass-weighting in amu, the
    // reduced masses computed will also be in amu.

    srs::size_t natoms  = static_cast<srs::size_t>(rot.atoms.size());
    srs::size_t natoms3 = 3 * natoms;
    dmat.resize(3, natoms, natoms3, 0.0);

    trans_rot_vec(dmat, n_tr_rot);

    // Construct n_vib other orthogonal vectors and put them at the beginning
    // of the transformation matrix D:

    srs::dmatrix tmp(natoms3, natoms3, dmat.data());
    srs::schmidt(tmp, n_tr_rot);
    dmat = srs::dcube(3, natoms, natoms3, tmp.data());
    shuffle(dmat, n_tr_rot);

    // Transform Cartesian Hessians to internal coordinates:

    tmp = srs::dmatrix(natoms3, natoms3, dmat.data());

    srs::dmatrix fc_int = get_mw_hessians() * tmp;

    lmat = srs::transpose(tmp) * fc_int;
}

void Molvib::print_cart_freq(std::ostream& to) const
{
    srs::Format<char> line;
    line.width(24).fill('-');
    srs::Format<double> fix;
    fix.fixed().width(10).precision(4);

    to << "Low frequencies (cm^-1):\n" << line('-') << '\n';

    auto cart_freqs  = calc_cart_freqs();
    srs::size_t n    = cart_freqs.size() > 9 ? 9 : cart_freqs.size();
    srs::size_t iter = 0;
    for (srs::size_t i = 0; i < n; ++i) {
        to << fix(cart_freqs(i)) << "  ";
        if (iter >= 5) {
            to << '\n';
            iter = 0;
        }
        ++iter;
    }
    to << '\n';
}

void Molvib::print_norm_modes(std::ostream& to) const
{
    srs::Format<char> line;
    line.fill('-').width(24);

    srs::Format<double> fix;
    fix.fixed().width(16).precision(4);

    to << "\nMode\tWavenumber/cm^-1\n" << line('-') << '\n';
    for (srs::size_t i = 0; i < freqs.size(); ++i) {
        to << i + 1 << '\t' << fix(freqs(i)) << '\n';
    }

    to << "\nZero-point vibrational energy: "
       << zero_point_energy() * datum::icm_to_kJ << " kJ/mol\n\n";

    line.width(21);
    to << "Normal mode analysis:\n" << line('-') << "\n\n";

    line.width(30);
    fix.width(10);

    srs::Format<double> fix2;
    fix2.fixed().width(6).precision(2);

    for (srs::size_t i = 0; i < freqs.size(); ++i) {
        to << "Mode:\t" << i + 1 << '\n'
           << line('-') << '\n'
           << "Wavenumber:\t" << fix(freqs(i)) << " cm^-1\n"
           << "Red. mass:\t" << fix(mu_freqs(i)) << " amu\n"
           << "Force const:\t" << fix(k_fc(i)) << " mDyne/A\n\n"
           << "Cartesian normal coordinates:\n"
           << line('-') << '\n'
           << "Atom\t  X\t  Y\t  Z\n";
        for (std::size_t j = 0; j < rot.atoms.size(); ++j) {
            to << rot.atoms[j].atomic_symbol;
            for (srs::size_t k = 0; k < 3; ++k) {
                to << '\t' << fix2(l_cart(k, j, i));
            }
            to << '\n';
        }
        to << '\n';
    }
}

void Molvib::shuffle(srs::dcube& dmat, srs::size_t n_tr_rot) const
{
    srs::dcube tmp(dmat);
    srs::size_t natoms3 = tmp.depths();
    srs::size_t n_vib   = natoms3 - n_tr_rot;

    for (srs::size_t k = 0; k < n_vib; ++k) {
        dmat.depth(k) = srs::dmatrix(tmp.depth(n_tr_rot + k));
    }
    for (srs::size_t k = n_vib; k < natoms3; ++k) {
        dmat.depth(k) = srs::dmatrix(tmp.depth(k - n_vib));
    }
}

void Molvib::freqs_unit_conv(srs::dvector& vib) const
{
    using namespace datum;

    const double factor  // conversion from atomic units to cm-1
        = 0.1 * (N_A * E_h / (4.0 * std::pow(pi * c_0 * a_0 * 1.0e-10, 2.0)));

    for (auto& vi : vib) {
        auto x = vi * factor;
        vi     = srs::sign(std::sqrt(std::abs(x)), x);
    }
}
