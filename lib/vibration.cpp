// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/vibration.h>
#include <numlib/constants.h>
#include <numlib/math.h>
#include <numlib/traits.h>
#include <cmath>

Chem::Vibration::Vibration(std::istream& from,
                           const std::string& key,
                           const std::vector<Chem::Element>& at,
                           const Numlib::Mat<double>& x,
                           const Numlib::Mat<double>& p)
    : atms(at), xyz(x), paxis(p)
{
    using namespace Stdutils;

    std::string token;
    Numlib::Vec<double> tmp;

    auto pos = find_token(from, key);
    if (pos != -1) {
        while (from >> token) {
            if (token == "frequencies") {
                from >> freqs;
                break;
            }
        }
    }
    pos = find_token(from, key);
    if (pos != -1) {
        while (from >> token) {
            if (token == "hessians") {
                from >> tmp;
                hess = Numlib::Symm_mat<double, Numlib::lo>(tmp);
                calc_normal_modes();
                break;
            }
        }
    }
}

double Chem::Vibration::zero_point_energy() const
{
    double zpe = 0.0;
    for (auto v : freqs) {
        if (v < 0.0) { // ignore imaginary frequencies
            continue;
        }
        else {
            zpe += v;
        }
    }
    return 0.5 * zpe;
}

void Chem::Vibration::calc_normal_modes()
{
    // Transform Cartesian force constants to internal coordinates:

    Numlib::Cube<double> dmat; // D matrix
    Numlib::Mat<double> lmat;  // L matrix
    int n_tr_rot;              // number of translational and rotational modes

    trans_hess_int_coord(dmat, lmat, n_tr_rot);
    Numlib::Symm_mat<double, Numlib::lo> fc_int(lmat);

    // Calculate vibrational frequencies:

    Index natoms = narrow_cast<Index>(atms.size());
    Index natoms3 = 3 * natoms;
    Index n_vib = natoms3 - n_tr_rot;

    freqs.resize(n_vib);
    lmat.resize(n_vib, n_vib);

    Numlib::eigs(fc_int, lmat, freqs);
    freqs_unit_conv(freqs);

    // Calculate reduced masses, force constants, and Cartesian displacements:

    using namespace Numlib::Constants;

    const double factor =
        4.0 * std::pow(c_0 * 100.0, 2.0) * m_u * pi * pi / 100.0;

    mu_freqs.resize(n_vib);
    k_fc.resize(n_vib);

    Numlib::Mat<double> tmp(3, natoms);
    Numlib::Mat<double> ltmp(natoms3, n_vib);

    for (Index i = 0; i < n_vib; ++i) {
        double y = 0.0;
        for (Index k1 = 0; k1 < natoms; ++k1) {
            for (int k2 = 0; k2 < 3; ++k2) {
                double x = 0.0;
                for (Index j = 0; j < n_vib; ++j) {
                    x += dmat(j, k1, k2) * lmat(j, i);
                }
                x /= std::sqrt(atms[k1].atomic_mass);
                y += x * x;
                tmp(k2, k1) = x;
            }
        }
        mu_freqs(i) = 1.0 / y;
        k_fc(i) = freqs(i) * freqs(i) * mu_freqs(i) * factor;
        y = 1.0 / std::sqrt(y);
        for (Index it = 0; it < natoms3; ++it) {
            ltmp(it, i) = y * tmp.data()[it];
        }
    }
    Numlib::Matrix_slice<3> ms(3, natoms, n_vib);
    l_cart = Numlib::Cube<double>(ms, ltmp.data());
}

void Chem::Vibration::print(std::ostream& to) const
{
    Stdutils::Format<char> line;
    line.width(26).fill('-');

    Stdutils::Format<double> fix;
    fix.fixed().width(8).precision(2);

    if (freqs.size() > 0) {
        Index it = 0;
        to << "Vibrational modes (cm^-1):\n" << line('-') << '\n';
        for (Index i = 0; i < freqs.size(); ++i) {
            to << fix(freqs(i));
            if ((it == 8) && (freqs.size() > 9)) {
                to << '\n';
            }
            it += 1;
        }
        double zpe = zero_point_energy();
        to << "\n\nZero-point vibrational energy: "
           << zpe / Numlib::Constants::au_to_icm << " Hartree\n\n";
    }
    if (!hess.empty()) {
        print_normal_modes(to);
    }
}

Numlib::Mat<double> Chem::Vibration::mw_hessians() const
{
    Numlib::Mat<double> hess_mw(hess.rows(), hess.cols());
    if (!hess.empty()) {
        for (Index i = 0; i < hess_mw.rows(); ++i) {
            for (Index j = 0; j < i + 1; ++j) {
                hess_mw(i, j) = hess(i, j);
                hess_mw(j, i) = hess_mw(i, j);
            }
        }
        for (Index i = 0; i < hess_mw.rows(); ++i) {
            for (Index j = 0; j < hess_mw.cols(); ++j) {
                hess_mw(i, j) /= std::sqrt(atms[i / 3].atomic_mass *
                                           atms[j / 3].atomic_mass);
            }
        }
    }
    return hess_mw;
}

void Chem::Vibration::trans_rot_vec(Numlib::Cube<double>& dmat,
                                    int& n_tr_rot) const
{
    // This function generates the vectors corresponding to translations and
    // infinitesimal rotations. The variable n_tr_rot is set to 5 or 6,
    // depending on how many such vectors there are. The function does not
    // take the symmetry of the molecule into account.

    Index natoms = narrow_cast<Index>(atms.size());

    // Determine center of mass, moments of inertia and rotation generators:

    auto xyz_ = xyz / Numlib::Constants::a_0;

    // Generate transformation matrix:

    for (Index i = 0; i < natoms; ++i) {
        double m = std::sqrt(atms[i].atomic_mass);
        double cx = xyz_(i, 0);
        double cy = xyz_(i, 1);
        double cz = xyz_(i, 2);
        double cxp = cx * paxis(0, 0) + cy * paxis(1, 0) + cz * paxis(2, 0);
        double cyp = cx * paxis(0, 1) + cy * paxis(1, 1) + cz * paxis(2, 1);
        double czp = cx * paxis(0, 2) + cy * paxis(1, 2) + cz * paxis(2, 2);
        dmat(0, i, 0) = m;
        dmat(1, i, 1) = m;
        dmat(2, i, 2) = m;
        dmat(3, i, 0) = (cyp * paxis(0, 2) - czp * paxis(0, 1)) * m;
        dmat(3, i, 1) = (cyp * paxis(1, 2) - czp * paxis(1, 1)) * m;
        dmat(3, i, 2) = (cyp * paxis(2, 2) - czp * paxis(2, 1)) * m;
        dmat(4, i, 0) = (czp * paxis(0, 0) - cxp * paxis(0, 2)) * m;
        dmat(4, i, 1) = (czp * paxis(1, 0) - cxp * paxis(1, 2)) * m;
        dmat(4, i, 2) = (czp * paxis(2, 0) - cxp * paxis(2, 2)) * m;
        dmat(5, i, 0) = (cxp * paxis(0, 1) - cyp * paxis(0, 0)) * m;
        dmat(5, i, 1) = (cxp * paxis(1, 1) - cyp * paxis(1, 0)) * m;
        dmat(5, i, 2) = (cxp * paxis(2, 1) - cyp * paxis(2, 0)) * m;
    }
    const double cutoff = 1.0e-12;

    n_tr_rot = 6;
    int i = 1;
    while (i < n_tr_rot) {
        auto ri = dmat.row(i - 1);
        double x = std::inner_product(ri.begin(), ri.end(), ri.begin(), 0.0);
        if (x < cutoff) {
            --n_tr_rot;
            ri = dmat.row(n_tr_rot);
        }
        else {
            ri *= (1.0 / std::sqrt(x));
            ++i;
        }
    }
    if (natoms == 1) {
        Assert::dynamic(n_tr_rot == 3, "bad n_tr_rot");
    }
}

void Chem::Vibration::trans_hess_int_coord(Numlib::Cube<double>& dmat,
                                           Numlib::Mat<double>& lmat,
                                           int& n_tr_rot) const
{
    // Set up coordinate vectors for translation and rotation about principal
    // axes of inertia.

    // Since these coordinate vectors involve mass-weighting in amu, the
    // reduced masses computed will also be in amu.

    Index natoms = narrow_cast<Index>(atms.size());
    Index natoms3 = 3 * natoms;

    dmat.resize(natoms3, natoms, 3);
    dmat = 0.0;

    trans_rot_vec(dmat, n_tr_rot);

    // Construct n_vib other orthogonal vectors and put them at the beginning
    // of the transformation matrix D:

    Numlib::Matrix_slice<2> ms2(natoms3, natoms3);
    Numlib::Mat<double> tmp(ms2, dmat.data());

    Numlib::schmidt(tmp, n_tr_rot);

    Numlib::Matrix_slice<3> ms3(natoms3, natoms, 3);
    dmat = Numlib::Cube<double>(ms3, tmp.data());

    shuffle(dmat, n_tr_rot);

    // Transform Cartesian Hessians to internal coordinates:

    tmp = Numlib::Mat<double>(ms2, dmat.data());

    auto fc_int = mw_hessians() * Numlib::transpose(tmp);
    lmat = tmp * fc_int;
}

void Chem::Vibration::shuffle(Numlib::Cube<double>& dmat, int n_tr_rot) const
{
    Numlib::Cube<double> tmp(dmat);

    Index natoms3 = tmp.rows();
    Index n_vib = natoms3 - n_tr_rot;

    for (Index k = 0; k < n_vib; ++k) {
        dmat.row(k) = tmp.row(n_tr_rot + k);
    }
    for (Index k = n_vib; k < natoms3; ++k) {
        dmat.row(k) = tmp.row(k - n_vib);
    }
}

void Chem::Vibration::freqs_unit_conv(Numlib::Vec<double>& vib) const
{
    using namespace Numlib::Constants;

    const double factor // conversion from atomic units to cm^-1
        = 0.1 * (N_A * E_h / (4.0 * std::pow(pi * c_0 * a_0 * 1.0e-10, 2.0)));

    for (auto& vi : vib) {
        auto x = vi * factor;
        vi = Numlib::sign(std::sqrt(std::abs(x)), x);
    }
}

void Chem::Vibration::print_cart_freqs(std::ostream& to) const
{
    Stdutils::Format<char> line;
    line.width(24).fill('-');
    Stdutils::Format<double> fix;
    fix.fixed().width(10).precision(4);

    auto h_mw = Numlib::Symm_mat<double, Numlib::lo>(mw_hessians());

    Numlib::Vec<double> hv_cart(h_mw.rows());
    Numlib::Mat<double> v;

    Numlib::eigs(h_mw, v, hv_cart);
    freqs_unit_conv(hv_cart);

    to << "Low frequencies (cm^-1):\n" << line('-') << '\n';

    Index n = hv_cart.size() > 9 ? 9 : hv_cart.size();
    Index iter = 0;
    for (Index i = 0; i < n; ++i) {
        to << fix(hv_cart(i)) << "  ";
        if (iter >= 5) {
            to << '\n';
            iter = 0;
        }
        ++iter;
    }
    to << '\n';
}

void Chem::Vibration::print_normal_modes(std::ostream& to) const
{
    print_cart_freqs(to);

    Stdutils::Format<char> line;
    line.fill('-').width(24);

    Stdutils::Format<double> fix;
    fix.fixed().width(16).precision(4);

    to << "\nMode\tWavenumber/cm^-1\n" << line('-') << '\n';
    for (Index i = 0; i < freqs.size(); ++i) {
        to << i + 1 << '\t' << fix(freqs(i)) << '\n';
    }

    to << "\nZero-point vibrational energy: "
       << zero_point_energy() * Numlib::Constants::icm_to_kJ << " kJ/mol\n\n";

    line.width(21);
    to << "Normal mode analysis:\n" << line('-') << "\n\n";

    line.width(30);
    fix.width(10);

    Stdutils::Format<double> fix2;
    fix2.fixed().width(6).precision(2);

    for (Index i = 0; i < freqs.size(); ++i) {
        to << "Mode:\t" << i + 1 << '\n'
           << line('-') << '\n'
           << "Wavenumber:\t" << fix(freqs(i)) << " cm^-1\n"
           << "Red. mass:\t" << fix(mu_freqs(i)) << " amu\n"
           << "Force const:\t" << fix(k_fc(i)) << " mDyne/A\n\n"
           << "Cartesian normal coordinates:\n"
           << line('-') << '\n'
           << "Atom\t  X\t  Y\t  Z\n";
        for (std::size_t j = 0; j < atms.size(); ++j) {
            to << atms[j].atomic_symbol;
            for (int k = 0; k < 3; ++k) {
                to << '\t' << fix2(l_cart(k, j, i));
            }
            to << '\n';
        }
        to << '\n';
    }
}

