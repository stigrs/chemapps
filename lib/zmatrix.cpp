// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/zmatrix.h>
#include <numlib/constants.h>
#include <numlib/math.h>
#include <numlib/traits.h>
#include <algorithm>
#include <cmath>

std::vector<Numlib::Vec<int>> Chem::Zmatrix::get_connectivities() const
{
    std::vector<Numlib::Vec<int>> connect(0);
    if (!atoms.empty()) {
        Numlib::Vec<int> ivec1 = {bond_connect(1)};
        connect.push_back(ivec1);
    }
    if (atoms.size() > 1) {
        Numlib::Vec<int> ivec2 = {bond_connect(2), angle_connect(2)};
        connect.push_back(ivec2);
    }
    if (atoms.size() > 2) {
        for (Index i = 3; i < bond_connect.size(); ++i) {
            Numlib::Vec<int> ivec3 = {bond_connect(i), angle_connect(i),
                                      dihedral_connect(i)};
            connect.push_back(ivec3);
        }
    }
    return connect;
}

void Chem::Zmatrix::rotate_moiety(const std::vector<int>& moiety, double value)
{
    if (atoms.size() > 3) {
        for (auto& mi : moiety) {
            double phi = get_dihedral(mi);
            set_dihedral(mi, phi + value);
        }
    }
}

void Chem::Zmatrix::build_zmat()
{
    Numlib::Mat<double> dist_mat;
    Numlib::pdist_matrix(dist_mat, xyz);

    for (int atom = 1; atom < narrow_cast<int>(atoms.size()); ++atom) {
        auto dmat_row = dist_mat.row(atom);
        auto dist = dmat_row(Numlib::slice{0, atom});
        bond_connect(atom) = find_nearest_atom(dist);
        distances(atom) = Numlib::min(dist);
        if (atom >= 2) {
            Numlib::Vec<int> iatms(3);
            iatms(0) = atom;
            iatms(1) = bond_connect(iatms(0));
            iatms(2) = bond_connect(iatms(1));
            if (iatms(2) == iatms(1)) {
                iatms(2) = find_new_connection(
                    iatms, bond_connect(Numlib::slice{0, atom}));
            }
            angle_connect(atom) = iatms(2);
            auto ai = xyz.row(iatms(0));
            auto aj = xyz.row(iatms(1));
            auto ak = xyz.row(iatms(2));
            angles(atom) = Numlib::angle(ai, aj, ak);
        }
        if (atom >= 3) {
            Numlib::Vec<int> iatms(4);
            iatms(0) = atom;
            iatms(1) = bond_connect(iatms(0));
            iatms(2) = angle_connect(iatms(0));
            iatms(3) = angle_connect(iatms(1));
            auto tmp = iatms(Numlib::slice{0, 3});
            auto it = std::find(tmp.begin(), tmp.end(), iatms(3));
            if (it != tmp.end()) {
                iatms(3) = find_new_connection(
                    iatms, bond_connect(Numlib::slice{0, atom}));
            }
            dihedral_connect(atom) = iatms(3);
            auto ai = xyz.row(iatms(0));
            auto aj = xyz.row(iatms(1));
            auto ak = xyz.row(iatms(2));
            auto al = xyz.row(iatms(3));
            dihedrals(atom) = Numlib::dihedral(ai, aj, ak, al);
        }
    }
}

void Chem::Zmatrix::build_xyz()
{
    xyz.resize(atoms.size(), 3);
    for (int atom = 0; atom < narrow_cast<int>(atoms.size()); ++atom) {
        xyz.row(atom) = calc_position(atom);
    }
}

int Chem::Zmatrix::find_nearest_atom(const Numlib::Vec<double>& dist) const
{
    double dist_min = Numlib::min(dist);
    int nearest_atom = -1;

    for (Index i = 0; i < dist.size(); ++i) {
        if (std::abs(dist(i) - dist_min) < 1.0e-12) {
            nearest_atom = narrow_cast<int>(i);
            break;
        }
    }
    return nearest_atom;
}

int Chem::Zmatrix::find_new_connection(
    const Numlib::Vec<int>& iatms, const Numlib::Vec<int>& connectivity) const
{
    int connection = 0;
    for (Index idx = 1; idx < connectivity.size(); ++idx) {
        // clang-format off
        if (std::find(iatms.begin(), iatms.end(), idx) == iatms.end() &&
            std::find(iatms.begin(), iatms.end(), connectivity(idx)) != iatms.end()) {
            connection = narrow_cast<int>(idx);
        }
        // clang-format off
    }
    return connection;
}

Numlib::Vec<double> Chem::Zmatrix::calc_position(int i) const
{
    Numlib::Vec<double> pos(3);
    double dst = 0.0;
    if (i > 1) {
        int j = bond_connect(i);
        int k = angle_connect(i);
        int l = dihedral_connect(i);
        Numlib::Vec<int> tmp = {i, j, k};
        if ((k == l) && (i > 0)) {  // prevent doubles
            l = find_new_connection(tmp, bond_connect(Numlib::slice{0, i}));
        }
        auto avec = xyz.row(j);
        auto bvec = xyz.row(k);
        dst = distances(i);
        double ang = Numlib::degtorad(angles(i));
        double tor;
        Numlib::Vec<double> cvec;
        if (i == 2) {  // third atom will be in the same plane as first two
            tor  = 90.0 * Numlib::Constants::pi / 180.0;
            cvec = {0.0, 1.0, 0.0};
        }
        else {  // fourth+ atoms require dihedral angle
            tor  = Numlib::degtorad(dihedrals(i));
            cvec = xyz.row(l);
        }
        auto v1 = avec - bvec;
        auto v2 = avec - cvec;
        auto n  = Numlib::cross(v1, v2);
        auto nn = Numlib::cross(v1, n);
        n /= Numlib::norm(n);
        nn /= Numlib::norm(nn);
        n *= -std::sin(tor);
        nn *= std::cos(tor);
        auto v3 = n + nn;
        v3 /= Numlib::norm(v3);
        v3 *= dst * std::sin(ang);
        v1 /= Numlib::norm(v1);
        v1 *= dst * std::cos(ang);
        pos = avec + v3 - v1;
    }
    else if (i == 1) {  // second atom dst away from origin along Z axis
        int j = bond_connect(i);
        dst = distances(i);
        pos = {xyz(j, 0) + dst, xyz(j, 1), xyz(j, 2)};
    }
    else if (i == 0) {  // first atom at the origin
        pos = 0.0;
    }
    return pos;
}

