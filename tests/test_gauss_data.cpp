// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/gauss_data.h>
#include <numlib/matrix.h>
#include <numlib/math.h>
#include <stdutils/stdutils.h>
#include <catch2/catch.hpp>
#include <fstream>
#include <cmath>

TEST_CASE("test_gauss_data")
{
    std::ifstream from;
    Stdutils::fopen(from, "test_gauss_data.inp");

    Numlib::Vec<double> ans = {
        -1.89865925E-04, -1.76751570E-16, 8.04584647E-01,  -3.41586331E-16,
        1.98174810E-14,  6.35148526E-01,  9.49329625E-05,  1.81796550E-17,
        4.34562970E-17,  -8.85968626E-05, -7.48866947E-17, -4.02292324E-01,
        2.16559539E-01,  3.81354171E-17,  4.39228651E-01,  1.08057743E-16,
        3.37433438E-01,  -3.17574263E-01, -2.67533069E-17, -2.76996488E-01,
        3.00228271E-01,  9.49329625E-05,  2.43493199E-16,  2.46627097E-16,
        -6.33609993E-06, -6.64102566E-18, -2.75427250E-17, -8.85968626E-05,
        2.39595307E-16,  -4.02292324E-01, -2.16559539E-01, -1.52988164E-17,
        -3.69363277E-02, -6.04369497E-02, -2.65825471E-16, 4.39228651E-01,
        1.87664060E-16,  -3.37433438E-01, -3.17574263E-01, 2.75352650E-17,
        6.04369497E-02,  1.73459917E-02,  -2.17458099E-16, 2.76996488E-01,
        3.00228271E-01};
    Numlib::Symm_mat<double, Numlib::lo> hess;

    Chem::Gauss_data gauss(from, Chem::fchk);
    gauss.get_hessians(hess);

    CHECK(hess.size() == ans.size());
    for (Index i = 0; i < ans.size(); ++i) {
        CHECK(std::abs(hess.data()[i] - ans(i)) < 1.0e-12);
    }
}

