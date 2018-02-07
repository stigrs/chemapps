//////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////

#include <chem/molecule_io.h>
#include <chem/zmatrix.h>
#include <srs/math.h>
#include <gsl/gsl>

Zmatrix::Zmatrix(const Zmatrix& zmat) : atoms(zmat.atoms), xyz(zmat.xyz)
{
    distances = zmat.distances;
    angles    = zmat.angles;
    dihedrals = zmat.dihedrals;

    bond_connect     = zmat.bond_connect;
    angle_connect    = zmat.angle_connect;
    dihedral_connect = zmat.dihedral_connect;
}

std::vector<srs::ivector> Zmatrix::get_connectivities() const
{
    std::vector<srs::ivector> connect(0);
    if (!atoms.empty()) {
        srs::ivector ivec1 = {bond_connect(1)};
        connect.push_back(ivec1);
    }
    if (atoms.size() > 1) {
        srs::ivector ivec2 = {bond_connect(2), angle_connect(2)};
        connect.push_back(ivec2);
    }
    if (atoms.size() > 2) {
        for (int i = 3; i < bond_connect.size(); ++i) {
            srs::ivector ivec3
                = {bond_connect(i), angle_connect(i), dihedral_connect(i)};
            connect.push_back(ivec3);
        }
    }
    return connect;
}
void Zmatrix::rotate_moiety(const std::vector<int>& moiety, double value)
{
    if (atoms.size() > 3) {
        for (std::size_t i = 0; i < moiety.size(); ++i) {
            double phi = get_dihedral(moiety[i]);
            set_dihedral(moiety[i], phi + value);
        }
    }
}

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
    srs::dmatrix dist_mat;
    chem::pdist_matrix(dist_mat, xyz);

    for (std::size_t atom = 1; atom < atoms.size(); ++atom) {
        srs::dvector dist  = dist_mat.row(atom).head(atom);
        bond_connect(atom) = find_nearest_atom(dist);
        distances(atom)    = dist.min();
        if (atom >= 2) {
            srs::ivector iatms(3);
            iatms(0) = atom;
            iatms(1) = bond_connect(iatms(0));
            iatms(2) = bond_connect(iatms(1));
            if (iatms(2) == iatms(1)) {
                iatms(2) = find_new_connection(iatms, bond_connect.head(atom));
            }
            angle_connect(atom) = iatms(2);
            srs::dvector ai     = xyz.row(iatms(0));
            srs::dvector aj     = xyz.row(iatms(1));
            srs::dvector ak     = xyz.row(iatms(2));
            angles(atom)        = chem::angle(ai, aj, ak);
        }
        if (atom >= 3) {
            srs::ivector iatms(4);
            iatms(0) = atom;
            iatms(1) = bond_connect(iatms(0));
            iatms(2) = angle_connect(iatms(0));
            iatms(3) = angle_connect(iatms(1));
            // Continue here:
            auto tmp = iatms.head(3);
            auto it std::find(tmp.begin(), tmp.end(), iatms(3));
            if (it != tmp.end()) {
                iatms(3) = find_new_connection(iatms, bond_connect.head(atom));
            }
            dihedral_connect(atom) = iatms(3);
            srs::dvector ai        = xyz.row(iatms(0));
            srs::dvector aj        = xyz.row(iatms(1));
            srs::dvector ak        = xyz.row(iatms(2));
            srs::dvector al        = xyz.row(iatms(3));
            dihedrals(atom)        = chem::dihedral(ai, aj, ak, al);
        }
    }
}

void Zmatrix::build_xyz()
{
    xyz = arma::zeros<arma::mat>(atoms.size(), 3);
    for (std::size_t atom = 0; atom < atoms.size(); ++atom) {
        xyz.row(atom) = calc_position(atom);
    }
}

int Zmatrix::find_nearest_atom(const arma::rowvec& dist) const
{
    double dist_min  = dist.min();
    int nearest_atom = -1;

    for (arma::uword i = 0; i < dist.size(); ++i) {
        if (chem::approx_equal(dist(i), dist_min)) {
            nearest_atom = gsl::narrow<int>(i);
            break;
        }
    }
    return nearest_atom;
}

int Zmatrix::find_new_connection(const arma::ivec& iatms,
                                 const arma::ivec& connectivity) const
{
    int connection = 0;
    for (arma::uword idx = 1; idx < connectivity.size(); ++idx) {
        if ((!arma::any(iatms == idx))
            && arma::any(iatms == connectivity(idx))) {
            connection = gsl::narrow<int>(idx);
        }
    }
    return connection;
}

arma::rowvec Zmatrix::calc_position(arma::sword i) const
{
    arma::rowvec pos(3);
    double dst = 0.0;
    if (i > 1) {
        arma::sword j  = bond_connect(i);
        arma::sword k  = angle_connect(i);
        arma::sword l  = dihedral_connect(i);
        arma::ivec tmp = {i, j, k};
        if ((k == l) && (i > 0)) {  // prevent doubles
            l = find_new_connection(tmp, bond_connect.head(i));
        }
        arma::rowvec avec = xyz.row(j);
        arma::rowvec bvec = xyz.row(k);
        dst               = distances(i);
        double ang        = chem::degtorad(angles(i));
        double tor;
        arma::rowvec cvec;
        if (i == 2) {  // third atom will be in the same plane as first two
            tor  = 90.0 * arma::datum::pi / 180.0;
            cvec = {0.0, 1.0, 0.0};
        }
        else {  // fourth+ atoms require dihedral angle
            tor  = chem::degtorad(dihedrals(i));
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
    else if (i == 1) {  // second atom dst away from origin along Z axis
        arma::sword j = bond_connect(i);
        dst           = distances(i);
        pos           = {xyz(j, 0) + dst, xyz(j, 1), xyz(j, 2)};
    }
    else if (i == 0) {  // first atom at the origin
        pos.zeros();
    }
    return pos;
}
