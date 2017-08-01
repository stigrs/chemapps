/**
   @file zmatrix.cpp
   
   This file is part of ChemApps - A C++ Chemistry Toolkit
   
   Copyright (C) 2016-2017  Stig Rune Sellevag
   
   ChemApps is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
 
   ChemApps is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <chem/zmatrix.h>
#include <chem/math.h>
#include <chem/molecule_io.h>


void Zmatrix::load(std::istream& from)
{
    chem::read_zmat_format(from, 
                           atoms,
                           distances,
                           angles,
                           dihedrals,
                           bond_connect,
                           angle_connect,
                           dihedral_connect);
    build_xyz();
}

void Zmatrix::print(std::ostream& to)
{
    chem::print_zmat_format(to, 
                            atoms,
                            distances,
                            angles,
                            dihedrals,
                            bond_connect,
                            angle_connect,
                            dihedral_connect);
}

void Zmatrix::build_zmat()
{
    arma::mat dist_mat;
    chem::pdist_matrix(dist_mat, xyz); 

    for (int atom = 1; atom < atoms.size(); ++atom) {
        arma::rowvec dist  = dist_mat.row(atom).head(atom);
        bond_connect(atom) = find_nearest_atom(dist);
        distances(atom)    = dist.min();
        if (atom >= 2) {
            arma::ivec iatms(3);
            iatms(0) = atom;
            iatms(1) = bond_connect(iatms(0));
            iatms(2) = bond_connect(iatms(1));
            if (iatms(2) == iatms(1)) {
                iatms(2) = find_new_connection(iatms, bond_connect.head(atom));
            }
            angle_connect(atom) = iatms(2);
            arma::vec ai = arma::conv_to<arma::vec>::from(xyz.row(iatms(0)));
            arma::vec aj = arma::conv_to<arma::vec>::from(xyz.row(iatms(1)));
            arma::vec ak = arma::conv_to<arma::vec>::from(xyz.row(iatms(2)));
            angles(atom) = chem::angle(ai, aj, ak);
        }
        if (atom >= 3) {
            arma::ivec iatms(4);
            iatms(0) = atom;
            iatms(1) = bond_connect(iatms(0));
            iatms(2) = angle_connect(iatms(0));
            iatms(3) = angle_connect(iatms(1));
            arma::uvec tmp = arma::find(iatms.head(3) == iatms(3));
            if (! tmp.empty()) {
                iatms(3) = find_new_connection(iatms, bond_connect.head(atom));
            }
            dihedral_connect(atom) = iatms(3);
            arma::vec ai = arma::conv_to<arma::vec>::from(xyz.row(iatms(0)));
            arma::vec aj = arma::conv_to<arma::vec>::from(xyz.row(iatms(1)));
            arma::vec ak = arma::conv_to<arma::vec>::from(xyz.row(iatms(2)));
            arma::vec al = arma::conv_to<arma::vec>::from(xyz.row(iatms(3)));
            dihedrals(atom) = chem::dihedral(ai, aj, ak, al);
        } 
    }
}

void Zmatrix::build_xyz()
{
    xyz = arma::zeros<arma::mat>(atoms.size(), 3);
    for (int atom = 0; atom < atoms.size(); ++atom) {
        xyz.row(atom) = calc_position(atom);
    }
}
        
int Zmatrix::find_nearest_atom(const arma::rowvec& dist) const
{
    double dist_min = dist.min();
    int nearest_atom = -1;

    for (int i = 0; i < dist.size(); ++i) {
        if (chem::approx_equal(dist(i), dist_min)) {
            nearest_atom = i;
            break;
        }
    }
    return nearest_atom;
}

int Zmatrix::find_new_connection(const arma::ivec& iatms,
                                 const arma::ivec& connectivity) const
{
    int connection = 0;
    for (int idx = 1; idx < connectivity.size(); ++idx) {
        if ((! arma::any(iatms == idx)) && 
            arma::any(iatms == connectivity(idx))) {
            connection = idx;
        }
    }
    return connection;
}

arma::rowvec Zmatrix::calc_position(int i) const
{
    arma::rowvec pos(3);
    double dst = 0.0;
    if (i > 1) {
        int j = bond_connect(i);
        int k = angle_connect(i);
        int l = dihedral_connect(i);
        arma::ivec tmp = {i, j, k};
        if ((k == l) && (i > 0)) { // prevent doubles
            l = find_new_connection(tmp, bond_connect.head(i));
        }
        arma::rowvec avec = xyz.row(j);
        arma::rowvec bvec = xyz.row(k);
        dst = distances(i);
        double ang = chem::degtorad(angles(i));
        double tor;
        arma::rowvec cvec;
        if (i == 2) { // third atom will be in the same plane as first two
            tor = 90.0 * arma::datum::pi / 180.0;
            cvec = {0.0, 1.0, 0.0};
        }
        else { // fourth+ atoms require dihedral angle
            tor = chem::degtorad(dihedrals(i));
            cvec = xyz.row(l);
        }
        arma::rowvec v1 = avec - bvec;
        arma::rowvec v2 = avec - cvec;
        arma::rowvec n  = arma::cross(v1, v2);
        arma::rowvec nn = arma::cross(v1, n);
        n /= arma::norm(n);
        nn /= arma::norm(nn);
        n *= -std::sin(tor);
        nn *= std::cos(tor);
        arma::rowvec v3 = n + nn;
        v3 /= arma::norm(v3);
        v3 *= dst * std::sin(ang);
        v1 /= arma::norm(v1);
        v1 *= dst * std::cos(ang);
        pos = avec + v3 - v1;
    }
    else if (i == 1) { // second atom dst away from origin along Z axis
        int j = bond_connect(i);
        dst = distances(i);
        pos = {xyz(j,0) + dst, xyz(j,1), xyz(j,2)};
    }
    else if (i == 0) { // first atom at the origin
        pos.zeros();
    }
    return pos;
}
