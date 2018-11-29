// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/geometry.h>
#include <numlib/constants.h>
#include <numlib/math.h>
#include <numlib/traits.h>
#include <algorithm>
#include <cmath>

Chem::Geometry::Geometry(std::istream& from, const std::string& key)
    : atms(), xyz()
{
    using namespace Stdutils;

    auto pos = find_token(from, key);
    if (pos != -1) {
        pos = find_token(from, "geometry", pos);
        if (pos != -1) {
            read_xyz_format(from, atms, xyz, info);
            if (!atms.empty()) {
                distances.resize(atms.size());
                angles.resize(atms.size());
                dihedrals.resize(atms.size());
                bond_connect.resize(atms.size());
                angle_connect.resize(atms.size());
                dihedral_connect.resize(atms.size());
                build_zmat();
            }
        }
    }
    if (atms.empty()) {
        pos = find_token(from, key);
        if (pos != -1) {
            from.ignore();
            from.clear();
            load_zmat(from);
        }
    }
}

std::vector<Numlib::Vec<Index>> Chem::Geometry::get_connectivities() const
{
    std::vector<Numlib::Vec<Index>> connect(0);
    if (!atms.empty()) {
        Numlib::Vec<Index> ivec1 = {bond_connect(1)};
        connect.push_back(ivec1);
    }
    if (atms.size() > 1) {
        Numlib::Vec<Index> ivec2 = {bond_connect(2), angle_connect(2)};
        connect.push_back(ivec2);
    }
    if (atms.size() > 2) {
        for (Index i = 3; i < bond_connect.size(); ++i) {
            Numlib::Vec<Index> ivec3 = {bond_connect(i), angle_connect(i),
                                        dihedral_connect(i)};
            connect.push_back(ivec3);
        }
    }
    return connect;
}

void Chem::Geometry::rotate_moiety(const std::vector<Index>& moiety,
                                   double value)
{
    if (atms.size() > 3) {
        for (auto& mi : moiety) {
            double phi = get_dihedral(mi);
            set_dihedral(mi, phi + value);
        }
    }
}

void Chem::Geometry::build_zmat()
{
    Numlib::Mat<double> dist_mat;
    Numlib::pdist_matrix(dist_mat, xyz);

    for (Index atom = 1; atom < narrow_cast<Index>(atms.size()); ++atom) {
        auto dmat_row = dist_mat.row(atom);
        auto dist = dmat_row(Numlib::slice{0, atom});
        bond_connect(atom) = find_nearest_atom(dist);
        distances(atom) = Numlib::min(dist);
        if (atom >= 2) {
            Numlib::Vec<Index> iatms(3);
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
            Numlib::Vec<Index> iatms(4);
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

void Chem::Geometry::build_xyz()
{
    xyz.resize(atms.size(), 3);
    for (Index atom = 0; atom < narrow_cast<Index>(atms.size()); ++atom) {
        xyz.row(atom) = calc_position(atom);
    }
}

Index Chem::Geometry::find_nearest_atom(const Numlib::Vec<double>& dist) const
{
    double dist_min = Numlib::min(dist);
    Index nearest_atom = -1;

    for (Index i = 0; i < dist.size(); ++i) {
        if (std::abs(dist(i) - dist_min) < 1.0e-12) {
            nearest_atom = i;
            break;
        }
    }
    return nearest_atom;
}

Index Chem::Geometry::find_new_connection(
    const Numlib::Vec<Index>& iatms,
    const Numlib::Vec<Index>& connectivity) const
{
    Index connection = 0;
    for (Index idx = 1; idx < connectivity.size(); ++idx) {
        // clang-format off
        if (std::find(iatms.begin(), iatms.end(), idx) == iatms.end() &&
            std::find(iatms.begin(), iatms.end(), connectivity(idx)) != iatms.end()) {
            connection = idx;
        }
        // clang-format off
    }
    return connection;
}

Numlib::Vec<double> Chem::Geometry::calc_position(Index i) const
{
    Numlib::Vec<double> pos(3);
    double dst = 0.0;
    if (i > 1) {
        Index j = bond_connect(i);
        Index k = angle_connect(i);
        Index l = dihedral_connect(i);
        Numlib::Vec<Index> tmp = {i, j, k};
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
        Index j = bond_connect(i);
        dst = distances(i);
        pos = {xyz(j, 0) + dst, xyz(j, 1), xyz(j, 2)};
    }
    else if (i == 0) {  // first atom at the origin
        pos = 0.0;
    }
    return pos;
}

