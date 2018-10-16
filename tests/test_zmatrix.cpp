// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/molecule.h>
#include <numlib/matrix.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <cmath>
#include <vector>

TEST_CASE("test_zmatrix")
{
    using namespace Chem;
    using namespace Numlib;
    using namespace Stdutils;

    std::ifstream from;
    fopen(from, "test_zmatrix.inp");

    Molecule mol(from);

    SECTION("num_atoms") { CHECK(mol.num_atoms() == 14); }
#if 0
    SECTION("bond_distances")
    {
        Vec<double> bonds = {1.54, 1.54, 1.54, 1.09, 1.09, 1.09, 1.09,
                             1.09, 1.09, 1.09, 1.09, 1.09, 1.09};

        for (int i = 1; i < mol.size(); ++i) {
            CHECK(std::abs(mol.get_zmat().get_distance(i) - bonds(i - 1)) <
                  1.0e-12);
        }
    }

    SECTION("bond_angles")
    {
        Vec<double> angles = {110.0, 110.0, 110.0, 110.0, 110.0, 110.0,
                              110.0, 110.0, 110.0, 110.0, 110.0, 110.0};

        for (int i = 2; i < mol.size(); ++i) {
            CHECK(std::abs(mol.get_zmat().get_angle(i) - angles(i - 2)) <
                  1.0e-12);
        }
    }

    SECTION("dihedral_angles")
    {
        Vec<double> dihedrals = {180.0, 0.0,   120.0, -120.0, 120.0, -120.0,
                                 60.0,  -60.0, 0.0,   120.0,  -120.0};

        for (int i = 3; i < mol.size(); ++i) {
            CHECK(std::abs(mol.get_zmat().get_dihedral(i) - dihedrals(i - 3)) <
                  1.0e-12);
        }
    }
    SECTION("rotate_moiety")
    {
        Mat<double> ans = {{0.00000000e+00, 0.00000000e+00, 0.00000000e+00},
                           {1.54000000e+00, 0.00000000e+00, 0.00000000e+00},
                           {2.06671102e+00, 8.86109501e-17, 1.44712664e+00},
                           {2.24685680e+00, 1.44712664e+00, 1.94207310e+00},
                           {-3.72801956e-01, 6.27181400e-17, 1.02426496e+00},
                           {-3.72801956e-01, 8.87039473e-01, -5.12132478e-01},
                           {-3.72801956e-01, -8.87039473e-01, -5.12132478e-01},
                           {1.91280196e+00, -8.87039473e-01, -5.12132478e-01},
                           {1.91280196e+00, 8.87039473e-01, -5.12132478e-01},
                           {3.02776125e+00, -5.12132478e-01, 1.49406052e+00},
                           {1.36067235e+00, -5.12132478e-01, 2.10083125e+00},
                           {1.96127393e+00, 2.14776513e+00, 1.15744062e+00},
                           {1.62151810e+00, 1.62228626e+00, 2.81749906e+00},
                           {3.28860700e+00, 1.62228626e+00, 2.21072833e+00}};

        std::vector<int> moiety = {3, 9, 10};
        mol.get_zmat().rotate_moiety(moiety, 90.0);
        auto res = mol.get_xyz();

        CHECK(same_extents(res, ans));

        for (Index i = 0; i < ans.rows(); ++i) {
            for (Index j = 0; j < ans.cols(); ++j) {
                CHECK(std::abs(res(i, j) - ans(i, j)) < 1.0e-8);
            }
        }
    }
#endif
}
