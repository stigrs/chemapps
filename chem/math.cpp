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

#include <chem/math.h>
#include <iostream>

double chem::dihedral(const arma::vec& a,
                      const arma::vec& b,
                      const arma::vec& c,
                      const arma::vec& d)
{
    arma::vec ab = arma::normalise(b - a);
    arma::vec bc = arma::normalise(c - b);
    arma::vec cd = arma::normalise(d - c);
    arma::vec n1 = arma::cross(ab, bc);
    arma::vec n2 = arma::cross(bc, cd);
    arma::vec m  = arma::cross(n1, bc);
    double x     = arma::dot(n1, n2);
    double y     = arma::dot(m, n2);

    double tau = radtodeg(std::atan2(y, x));
    if (std::abs(tau) < 1.0e-8) {  // avoid very small angles close to zero
        tau = 0.0;
    }
    return tau;
}

void chem::pdist_matrix(arma::mat& dm, const arma::mat& mat)
{
    dm = arma::zeros<arma::mat>(mat.n_rows, mat.n_rows);

    arma::rowvec dij(3);

    for (arma::uword j = 0; j < dm.n_cols; ++j) {
        for (arma::uword i = j; i < dm.n_rows; ++i) {
            if (i != j) {
                dij = mat.row(i) - mat.row(j);
                dm(i, j) = arma::norm(dij);
                dm(j, i) = dm(i, j);
            }
        }
    }
}

void chem::translate(arma::mat& xyz, double dx, double dy, double dz)
{
    for (arma::uword i = 0; i < xyz.n_rows; ++i) {
        xyz(i, 0) += dx;
        xyz(i, 1) += dy;
        xyz(i, 2) += dz;
    }
}

void chem::rotate(arma::mat& xyz, const arma::mat33& rotm)
{
    arma::vec3 xold;
    arma::vec3 xnew;

    for (arma::uword i = 0; i < xyz.n_rows; ++i) {
        arma::vec3 xyz_new = rotm * xyz.row(i).t();
        xyz(i, 0) = xyz_new(0);
        xyz(i, 1) = xyz_new(1);
        xyz(i, 2) = xyz_new(2);
    }
}
