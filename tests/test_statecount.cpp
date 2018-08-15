////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2018 Stig Rune Sellevag. All rights reserved.
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

#include <chem/energy_levels.h>
#include <chem/statecount.h>
#include <srs/array.h>
#include <srs/datum.h>
#include <srs/math.h>
#include <catch/catch.hpp>


TEST_CASE("test_statecount")
{
    namespace sc = statecount;

    SECTION("Table_4.1-4.2")  // Holbrook, Pilling and Robertson (1996)
    {
        srs::dvector vibr = {1500.0, 1200.0, 600.0};
        srs::dvector wans
            = {1.0, 1.0, 2.0, 2.0, 4.0, 5.0, 7.0, 8.0, 11.0, 13.0, 17.0};
        srs::dvector dans = {3.33e-3,
                             0.0,
                             3.33e-3,
                             0.0,
                             6.67e-3,
                             3.33e-3,
                             6.67e-3,
                             2.22e-3,
                             1.0e-2,
                             6.67e-3,
                             1.33e-2};

        double emin   = 0.0;
        double emax   = 3000.0;
        double egrain = 300.0;
        int ngrains   = 1 + srs::round<int>((emax - emin) / egrain);
        auto wres     = sc::bswine(vibr, ngrains, egrain, true);
        auto dres     = sc::bswine(vibr, ngrains, egrain);

        CHECK(srs::approx_equal(wres, wans, 1.0e-12));
        CHECK(srs::approx_equal(dres, dans, 1.0e-2, "reldiff"));
    }

    SECTION("Ethane_free_rotor")  // Stein and Rabinovitch (1973)
    {
        double emax   = 40000.0;
        double egrain = 5.0;
        int ngrains   = 1 + srs::round<int>(emax / egrain);

        srs::dvector vibr = {2915.0,
                             2915.0,
                             1388.0,
                             995.0,
                             1370.0,
                             2974.0,
                             2974.0,
                             1460.0,
                             1460.0,
                             822.0,
                             822.0,
                             2950.0,
                             2950.0,
                             1469.0,
                             1469.0,
                             1190.0,
                             1190.0};

        double sigma = 3.0;
        double rotc  = 10.704;

        auto wrot = sc::free_rotor(sigma, rotc, ngrains, egrain, true);
        auto wsum = sc::bswine(vibr, ngrains, egrain, true, wrot);

        srs::dvector en = {500.0,
                           1000.0,
                           2000.0,
                           4000.0,
                           6000.0,
                           8000.0,
                           10000.0,
                           20000.0,
                           30000.0,
                           40000.0};

        srs::dvector wsum_ans = {4.3,
                                 12.7,
                                 88.0,
                                 1728.0,
                                 1.8211e+4,
                                 13.071e+4,
                                 7.1416e+5,
                                 5.0222e+8,
                                 5.0613e+10,
                                 1.86446e+12};

        srs::ivector idx;
        for (auto ei : en) {
            idx.push_back(srs::round<int>(ei / egrain));
        }
        for (int i = 0; i < idx.size(); ++i) {
            CHECK(srs::approx_equal(
                wsum(idx(i)), wsum_ans(i), 6.0e-2, "reldiff"));
        }
    }

    SECTION("Ethane_hindered_rotor")  // Stein and Rabinovitch (1973)
    {
        double emax   = 40000.0;
        double egrain = 5.0;
        int ngrains   = 1 + srs::round<int>(emax / egrain);

        srs::dvector vibr = {2915.0,
                             2915.0,
                             1388.0,
                             995.0,
                             1370.0,
                             2974.0,
                             2974.0,
                             1460.0,
                             1460.0,
                             822.0,
                             822.0,
                             2950.0,
                             2950.0,
                             1469.0,
                             1469.0,
                             1190.0,
                             1190.0};

        double sigma = 3.0;
        double rotc  = 10.704;
        double v0    = 1024.0;

        auto wrot = sc::hindered_rotor(sigma, rotc, v0, ngrains, egrain, true);
        auto wsum = sc::bswine(vibr, ngrains, egrain, true, wrot);

        srs::dvector en = {500.0,
                           1000.0,
                           2000.0,
                           4000.0,
                           6000.0,
                           8000.0,
                           10000.0,
                           20000.0,
                           30000.0,
                           40000.0};

        srs::dvector wsum_ans = {2.0,
                                 8.0,
                                 55.7,
                                 1187.0,
                                 1.3111e+4,
                                 9.8719e+4,
                                 5.7428e+5,
                                 4.2130e+8,
                                 4.4204e+10,
                                 1.66923e+12};

        srs::ivector idx;
        for (auto ei : en) {
            idx.push_back(srs::round<int>(ei / egrain));
        }
        for (int i = 0; i < idx.size(); ++i) {
            CHECK(
                srs::approx_equal(wsum(idx(i)), wsum_ans(i), 0.23, "reldiff"));
        }
    }

    SECTION("NH3_bswine")  // Stein and Rabinovitch (1973)
    {
        srs::dvector vibr = {3337.0, 950.0, 3414.0, 3414.0, 1628.0, 1628.0};

        double emax   = 34976.0;
        double egrain = 1.0;
        int ngrains   = 1 + srs::round<int>(emax / egrain);

        auto wsum = sc::bswine(vibr, ngrains, egrain, true);

        srs::dvector en = {3.4980E+03,
                           6.9950E+03,
                           1.0493E+04,
                           1.399E+04,
                           1.7488E+04,
                           2.4483E+04};

        srs::dvector wsum_ans = {14., 94., 375., 1135., 2916., 13518.};

        srs::ivector idx;
        for (auto ei : en) {
            idx.push_back(srs::round<int>(ei / egrain));
        }
        for (int i = 0; i < idx.size(); ++i) {
            CHECK(srs::approx_equal(
                wsum(idx(i)), wsum_ans(i), 1.0e-8, "reldiff"));
        }
    }

    SECTION("Cyclopropane_bswine")  // Stein and Rabinovitch (1973)
    {
        srs::dvector vibr = {3221., 3221., 3221., 3221., 3221., 3221., 1478.,
                             1478., 1478., 1118., 1118., 1118., 1118., 1118.,
                             1118., 1118., 879.,  879.,  879.,  750.,  750.};

        double emax   = 34976.0;
        double egrain = 1.0;
        int ngrains   = 1 + srs::round<int>(emax / egrain);

        auto wsum = sc::bswine(vibr, ngrains, egrain, true);

        srs::dvector en = {3.4980E+03,
                           6.9950E+03,
                           1.0493E+04,
                           1.399E+04,
                           1.7488E+04,
                           2.4483E+04,
                           3.4976E+04};

        srs::dvector wsum_ans = {802.,
                                 77522.,
                                 2680083.,
                                 49612574.,
                                 6.11428E+08,
                                 4.0751802E+10,
                                 5.8325046E+12};

        srs::ivector idx;
        for (auto ei : en) {
            idx.push_back(srs::round<int>(ei / egrain));
        }
        for (int i = 0; i < idx.size(); ++i) {
            CHECK(srs::approx_equal(
                wsum(idx(i)), wsum_ans(i), 1.0e-8, "reldiff"));
        }
    }

    SECTION("NH3_steinrab")  // Stein and Rabinovitch (1973)
    {
        srs::dvector vibr = {3337.0, 950.0, 3414.0, 3414.0, 1628.0, 1628.0};

        double emax   = 34976.0;
        double egrain = 1.0;
        int ngrains   = 1 + srs::round<int>(emax / egrain);

        auto wsum = sc::steinrab(vibr, 0.0, 0.0, 0.0, ngrains, egrain, true);

        srs::dvector en = {3.4980E+03,
                           6.9950E+03,
                           1.0493E+04,
                           1.399E+04,
                           1.7488E+04,
                           2.4483E+04};

        srs::dvector wsum_ans = {14., 94., 375., 1135., 2916., 13518.};

        srs::ivector idx;
        for (auto ei : en) {
            idx.push_back(srs::round<int>(ei / egrain));
        }
        for (int i = 0; i < idx.size(); ++i) {
            CHECK(srs::approx_equal(
                wsum(idx(i)), wsum_ans(i), 1.0e-8, "reldiff"));
        }
    }

    SECTION("Cyclopropane_steinrab")  // Stein and Rabinovitch (1973)
    {
        srs::dvector vibr = {3221., 3221., 3221., 3221., 3221., 3221., 1478.,
                             1478., 1478., 1118., 1118., 1118., 1118., 1118.,
                             1118., 1118., 879.,  879.,  879.,  750.,  750.};

        double emax   = 34976.0;
        double egrain = 1.0;
        int ngrains   = 1 + srs::round<int>(emax / egrain);

        auto wsum = sc::steinrab(vibr, 0.0, 0.0, 0.0, ngrains, egrain, true);

        srs::dvector en = {3.4980E+03,
                           6.9950E+03,
                           1.0493E+04,
                           1.399E+04,
                           1.7488E+04,
                           2.4483E+04,
                           3.4976E+04};

        srs::dvector wsum_ans = {802.,
                                 77522.,
                                 2680083.,
                                 49612574.,
                                 6.11428E+08,
                                 4.0751802E+10,
                                 5.8325046E+12};

        srs::ivector idx;
        for (auto ei : en) {
            idx.push_back(srs::round<int>(ei / egrain));
        }
        for (int i = 0; i < idx.size(); ++i) {
            CHECK(srs::approx_equal(
                wsum(idx(i)), wsum_ans(i), 1.0e-8, "reldiff"));
        }
    }

    SECTION("Ethane_free_rotor_steinrab")  // Stein and Rabinovitch (1973)
    {
        double emax   = 40000.0;
        double egrain = 5.0;
        int ngrains   = 1 + srs::round<int>(emax / egrain);

        srs::dvector vibr = {2915.0,
                             2915.0,
                             1388.0,
                             995.0,
                             1370.0,
                             2974.0,
                             2974.0,
                             1460.0,
                             1460.0,
                             822.0,
                             822.0,
                             2950.0,
                             2950.0,
                             1469.0,
                             1469.0,
                             1190.0,
                             1190.0};

        double sigma = 3.0;
        double rotc  = 10.704;

        auto wsum = sc::steinrab(vibr, sigma, rotc, 0.0, ngrains, egrain, true);

        srs::dvector en = {500.0,
                           1000.0,
                           2000.0,
                           4000.0,
                           6000.0,
                           8000.0,
                           10000.0,
                           20000.0,
                           30000.0,
                           40000.0};

        srs::dvector wsum_ans = {4.3,
                                 12.7,
                                 88.0,
                                 1728.0,
                                 1.8211e+4,
                                 13.071e+4,
                                 7.1416e+5,
                                 5.0222e+8,
                                 5.0613e+10,
                                 1.86446e+12};

        srs::ivector idx;
        for (auto ei : en) {
            idx.push_back(srs::round<int>(ei / egrain));
        }
        for (int i = 0; i < idx.size(); ++i) {
            CHECK(srs::approx_equal(
                wsum(idx(i)), wsum_ans(i), 5.0e-2, "reldiff"));
        }
    }

    SECTION("Ethane_hindered_rotor_steinrab")  // Stein and Rabinovitch (1973)
    {
        double emax   = 40000.0;
        double egrain = 5.0;
        int ngrains   = 1 + srs::round<int>(emax / egrain);

        srs::dvector vibr = {2915.0,
                             2915.0,
                             1388.0,
                             995.0,
                             1370.0,
                             2974.0,
                             2974.0,
                             1460.0,
                             1460.0,
                             822.0,
                             822.0,
                             2950.0,
                             2950.0,
                             1469.0,
                             1469.0,
                             1190.0,
                             1190.0};

        double sigma = 3.0;
        double rotc  = 10.704;
        double v0    = 1024.0;

        auto wsum = sc::steinrab(vibr, sigma, rotc, v0, ngrains, egrain, true);

        srs::dvector en = {500.0,
                           1000.0,
                           2000.0,
                           4000.0,
                           6000.0,
                           8000.0,
                           10000.0,
                           20000.0,
                           30000.0,
                           40000.0};

        srs::dvector wsum_ans = {2.0,
                                 8.0,
                                 55.7,
                                 1187.0,
                                 1.3111e+4,
                                 9.8719e+4,
                                 5.7428e+5,
                                 4.2130e+8,
                                 4.4204e+10,
                                 1.66923e+12};

        srs::ivector idx;
        for (auto ei : en) {
            idx.push_back(srs::round<int>(ei / egrain));
        }
        for (int i = 0; i < idx.size(); ++i) {
            CHECK(
                srs::approx_equal(wsum(idx(i)), wsum_ans(i), 0.17, "reldiff"));
        }
    }
}
