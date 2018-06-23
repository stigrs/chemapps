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

#include <chem/molecule.h>
#include <srs/datum.h>
#include <srs/packed.h>
#include <srs/utils.h>
#include <catch/catch.hpp>
#include <cmath>
#include <fstream>
#include <iostream>


TEST_CASE("test_molvib")
{
    SECTION("H2O")
    {
        std::ifstream from;
        srs::fopen(from, "test_molvib_h2o.inp");

        Molecule mol(from);

        double zpe_ans = 0.024386;
        double zpe     = mol.get_vib().zero_point_energy() / datum::au_to_icm;
        CHECK(std::abs(zpe - zpe_ans) < 1.0e-6);

        srs::dvector nu_ans = {2169.7566, 4141.6012, 4392.7809};
        srs::dvector mu_ans = {1.0785, 1.0491, 1.0774};
        srs::dvector k_ans  = {2.9915, 10.6023, 12.2488};

        srs::dvector nu = mol.get_vib().get_freqs();
        srs::dvector mu = mol.get_vib().get_red_mass();
        srs::dvector k  = mol.get_vib().get_force_constant();

        CHECK(srs::approx_equal(nu, nu_ans, 1.0e-4));
        CHECK(srs::approx_equal(mu, mu_ans, 1.0e-4));
        CHECK(srs::approx_equal(k, k_ans, 1.0e-4));
    }

    SECTION("CO2")
    {
        std::ifstream from;
        srs::fopen(from, "test_molvib_co2.inp");

        Molecule mol(from);

        double zpe_ans = 0.011626;
        double zpe     = mol.get_vib().zero_point_energy() / datum::au_to_icm;
        CHECK(std::abs(zpe - zpe_ans) < 1.0e-6);

        srs::dvector nu_ans = {566.2442, 566.2442, 1435.2233, 2535.7089};
        srs::dvector mu_ans = {12.8774, 12.8774, 15.9949, 12.8774};
        srs::dvector k_ans  = {2.4327, 2.4327, 19.4120, 48.7838};

        srs::dvector nu = mol.get_vib().get_freqs();
        srs::dvector mu = mol.get_vib().get_red_mass();
        srs::dvector k  = mol.get_vib().get_force_constant();

        CHECK(srs::approx_equal(nu, nu_ans, 1.0e-4));
        CHECK(srs::approx_equal(mu, mu_ans, 1.0e-4));
        CHECK(srs::approx_equal(k, k_ans, 1.0e-4));
    }

    SECTION("CH4OH")
    {
        std::ifstream from;
        srs::fopen(from, "test_molvib_ch4oh.inp");

        Molecule mol(from);

        double zpe_ans = 0.051319;
        double zpe     = mol.get_vib().zero_point_energy() / datum::au_to_icm;
        CHECK(std::abs(zpe - zpe_ans) < 1.0e-6);

        srs::dvector nu_ans = {-1752.6530,
                               62.301,
                               326.8547,
                               332.5309,
                               775.5512,
                               909.5669,
                               1212.3977,
                               1275.3228,
                               1340.0812,
                               1451.6446,
                               1476.6141,
                               3103.4053,
                               3244.1604,
                               3247.8661,
                               3767.9277};
        srs::dvector mu_ans = {1.1201,
                               1.0445,
                               1.0905,
                               1.0684,
                               2.1719,
                               1.5377,
                               1.1215,
                               1.0879,
                               1.1208,
                               1.0421,
                               1.0216,
                               1.0242,
                               1.1068,
                               1.1074,
                               1.0669};
        srs::dvector k_ans  = {2.0272,
                              0.0024,
                              0.0686,
                              0.0696,
                              0.7697,
                              0.7495,
                              0.9713,
                              1.0425,
                              1.1859,
                              1.2939,
                              1.3124,
                              5.8116,
                              6.8635,
                              6.8826,
                              8.9245};

        srs::dvector nu = mol.get_vib().get_freqs();
        srs::dvector mu = mol.get_vib().get_red_mass();
        srs::dvector k  = mol.get_vib().get_force_constant();

        CHECK(srs::approx_equal(nu, nu_ans, 1.0e-4));
        CHECK(srs::approx_equal(mu, mu_ans, 1.0e-4));
        CHECK(srs::approx_equal(k, k_ans, 1.0e-4));
    }
}
