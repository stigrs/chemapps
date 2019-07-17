// Copyright (c) 2017 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#include <chem/energy_levels.h>
#include <chem/statecount.h>
#include <numlib/matrix.h>
#include <numlib/constants.h>
#include <numlib/math.h>
#include <catch2/catch.hpp>
#include <cmath>
#include <vector>

TEST_CASE("test_statecount")
{
    namespace Sc = Chem::Statecount;

    SECTION("Table_4.1-4.2") // Holbrook, Pilling and Robertson (1996)
    {
        Numlib::Vec<double> vibr = {1500.0, 1200.0, 600.0};
        Numlib::Vec<double> wans = {1.0, 1.0, 2.0,  2.0,  4.0, 5.0,
                                    7.0, 8.0, 11.0, 13.0, 17.0};
        Numlib::Vec<double> dans = {3.33e-3, 0.0,     3.33e-3, 0.0,
                                    6.67e-3, 3.33e-3, 6.67e-3, 3.33e-3,
                                    1.0e-2,  6.67e-3, 1.33e-2};

        double emin = 0.0;
        double emax = 3000.0;
        double egrain = 300.0;
        int ngrains = 1 + Numlib::round<int>((emax - emin) / egrain);
        auto wres = Sc::bswine(vibr, ngrains, egrain, true);
        auto dres = Sc::bswine(vibr, ngrains, egrain);

        for (Index i = 0; i < wans.size(); ++i) {
            CHECK(std::abs(wres(i) - wans(i)) < 1.0e-12);
        }
        for (Index i = 0; i < dans.size(); ++i) {
            CHECK(std::abs(dres(i) - dans(i)) < 5.0e-5);
        }
    }

    SECTION("Ethane_free_rotor") // Stein and Rabinovitch (1973)
    {
        double emax = 40000.0;
        double egrain = 5.0;
        int ngrains = 1 + Numlib::round<int>(emax / egrain);

        Numlib::Vec<double> vibr = {2915.0, 2915.0, 1388.0, 995.0,  1370.0,
                                    2974.0, 2974.0, 1460.0, 1460.0, 822.0,
                                    822.0,  2950.0, 2950.0, 1469.0, 1469.0,
                                    1190.0, 1190.0};

        double sigma = 3.0;
        double rotc = 10.704;

        auto wrot = Sc::free_rotor(sigma, rotc, ngrains, egrain, true);
        auto wsum = Sc::bswine(vibr, ngrains, egrain, true, wrot);

        Numlib::Vec<double> en = {500.0,  1000.0,  2000.0,  4000.0,  6000.0,
                                  8000.0, 10000.0, 20000.0, 30000.0, 40000.0};

        Numlib::Vec<double> wsum_ans = {
            4.3,       12.7,      88.0,      1728.0,     1.8211e+4,
            13.071e+4, 7.1416e+5, 5.0222e+8, 5.0613e+10, 1.86446e+12};

        std::vector<int> idx;
        for (auto ei : en) {
            idx.push_back(Numlib::round<int>(ei / egrain));
        }
        for (std::size_t i = 0; i < idx.size(); ++i) {
            CHECK(std::abs(wsum(idx[i]) - wsum_ans(i)) / wsum_ans(i) < 6.0e-2);
        }
    }

    SECTION("Ethane_hindered_rotor") // Stein and Rabinovitch (1973)
    {
        double emax = 40000.0;
        double egrain = 5.0;
        int ngrains = 1 + Numlib::round<int>(emax / egrain);

        Numlib::Vec<double> vibr = {2915.0, 2915.0, 1388.0, 995.0,  1370.0,
                                    2974.0, 2974.0, 1460.0, 1460.0, 822.0,
                                    822.0,  2950.0, 2950.0, 1469.0, 1469.0,
                                    1190.0, 1190.0};

        double sigma = 3.0;
        double rotc = 10.704;
        double v0 = 1024.0;

        auto wrot = Sc::hindered_rotor(sigma, rotc, v0, ngrains, egrain, true);
        auto wsum = Sc::bswine(vibr, ngrains, egrain, true, wrot);

        Numlib::Vec<double> en = {500.0,  1000.0,  2000.0,  4000.0,  6000.0,
                                  8000.0, 10000.0, 20000.0, 30000.0, 40000.0};

        Numlib::Vec<double> wsum_ans = {
            2.0,       8.0,       55.7,      1187.0,     1.3111e+4,
            9.8719e+4, 5.7428e+5, 4.2130e+8, 4.4204e+10, 1.66923e+12};

        std::vector<int> idx;
        for (auto ei : en) {
            idx.push_back(Numlib::round<int>(ei / egrain));
        }
        for (std::size_t i = 0; i < idx.size(); ++i) {
            CHECK(std::abs(wsum(idx[i]) - wsum_ans(i)) / wsum_ans(i) < 0.28);
        }
    }

    SECTION("NH3_bswine") // Stein and Rabinovitch (1973)
    {
        Numlib::Vec<double> vibr = {3337.0, 950.0,  3414.0,
                                    3414.0, 1628.0, 1628.0};

        double emax = 34976.0;
        double egrain = 1.0;
        int ngrains = 1 + Numlib::round<int>(emax / egrain);

        auto wsum = Sc::bswine(vibr, ngrains, egrain, true);

        Numlib::Vec<double> en = {3.4980E+03, 6.9950E+03, 1.0493E+04,
                                  1.399E+04,  1.7488E+04, 2.4483E+04};

        Numlib::Vec<double> wsum_ans = {14., 94., 375., 1135., 2916., 13518.};

        std::vector<int> idx;
        for (auto ei : en) {
            idx.push_back(Numlib::round<int>(ei / egrain));
        }
        for (std::size_t i = 0; i < idx.size(); ++i) {
            CHECK(std::abs(wsum(idx[i]) - wsum_ans(i)) / wsum_ans(i) < 1.0e-8);
        }
    }

    SECTION("Cyclopropane_bswine") // Stein and Rabinovitch (1973)
    {
        Numlib::Vec<double> vibr = {3221., 3221., 3221., 3221., 3221., 3221.,
                                    1478., 1478., 1478., 1118., 1118., 1118.,
                                    1118., 1118., 1118., 1118., 879.,  879.,
                                    879.,  750.,  750.};

        double emax = 34976.0;
        double egrain = 1.0;
        int ngrains = 1 + Numlib::round<int>(emax / egrain);

        auto wsum = Sc::bswine(vibr, ngrains, egrain, true);

        Numlib::Vec<double> en = {3.4980E+03, 6.9950E+03, 1.0493E+04, 1.399E+04,
                                  1.7488E+04, 2.4483E+04, 3.4976E+04};

        Numlib::Vec<double> wsum_ans = {
            802.,        77522.,        2680083.,     49612574.,
            6.11428E+08, 4.0751802E+10, 5.8325046E+12};

        std::vector<int> idx;
        for (auto ei : en) {
            idx.push_back(Numlib::round<int>(ei / egrain));
        }
        for (std::size_t i = 0; i < idx.size(); ++i) {
            CHECK(std::abs(wsum(idx[i]) - wsum_ans(i)) / wsum_ans(i) < 1.0e-8);
        }
    }

    SECTION("NH3_steinrab") // Stein and Rabinovitch (1973)
    {
        Numlib::Vec<double> vibr = {3337.0, 950.0,  3414.0,
                                    3414.0, 1628.0, 1628.0};

        double emax = 34976.0;
        double egrain = 1.0;
        int ngrains = 1 + Numlib::round<int>(emax / egrain);

        auto wsum = Sc::steinrab(vibr, 0.0, 0.0, 0.0, ngrains, egrain, true);

        Numlib::Vec<double> en = {3.4980E+03, 6.9950E+03, 1.0493E+04,
                                  1.399E+04,  1.7488E+04, 2.4483E+04};

        Numlib::Vec<double> wsum_ans = {14., 94., 375., 1135., 2916., 13518.};

        std::vector<int> idx;
        for (auto ei : en) {
            idx.push_back(Numlib::round<int>(ei / egrain));
        }
        for (std::size_t i = 0; i < idx.size(); ++i) {
            CHECK(std::abs(wsum(idx[i]) - wsum_ans(i)) / wsum_ans(i) < 1.0e-8);
        }
    }

    SECTION("Cyclopropane_steinrab") // Stein and Rabinovitch (1973)
    {
        Numlib::Vec<double> vibr = {3221., 3221., 3221., 3221., 3221., 3221.,
                                    1478., 1478., 1478., 1118., 1118., 1118.,
                                    1118., 1118., 1118., 1118., 879.,  879.,
                                    879.,  750.,  750.};

        double emax = 34976.0;
        double egrain = 1.0;
        int ngrains = 1 + Numlib::round<int>(emax / egrain);

        auto wsum = Sc::steinrab(vibr, 0.0, 0.0, 0.0, ngrains, egrain, true);

        Numlib::Vec<double> en = {3.4980E+03, 6.9950E+03, 1.0493E+04, 1.399E+04,
                                  1.7488E+04, 2.4483E+04, 3.4976E+04};

        Numlib::Vec<double> wsum_ans = {
            802.,        77522.,        2680083.,     49612574.,
            6.11428E+08, 4.0751802E+10, 5.8325046E+12};

        std::vector<int> idx;
        for (auto ei : en) {
            idx.push_back(Numlib::round<int>(ei / egrain));
        }
        for (std::size_t i = 0; i < idx.size(); ++i) {
            CHECK(std::abs(wsum(idx[i]) - wsum_ans(i)) / wsum_ans(i) < 1.0e-8);
        }
    }

    SECTION("Ethane_free_rotor_steinrab") // Stein and Rabinovitch (1973)
    {
        double emax = 40000.0;
        double egrain = 5.0;
        int ngrains = 1 + Numlib::round<int>(emax / egrain);

        Numlib::Vec<double> vibr = {2915.0, 2915.0, 1388.0, 995.0,  1370.0,
                                    2974.0, 2974.0, 1460.0, 1460.0, 822.0,
                                    822.0,  2950.0, 2950.0, 1469.0, 1469.0,
                                    1190.0, 1190.0};

        double sigma = 3.0;
        double rotc = 10.704;

        auto wsum = Sc::steinrab(vibr, sigma, rotc, 0.0, ngrains, egrain, true);

        Numlib::Vec<double> en = {500.0,  1000.0,  2000.0,  4000.0,  6000.0,
                                  8000.0, 10000.0, 20000.0, 30000.0, 40000.0};

        Numlib::Vec<double> wsum_ans = {
            4.3,       12.7,      88.0,      1728.0,     1.8211e+4,
            13.071e+4, 7.1416e+5, 5.0222e+8, 5.0613e+10, 1.86446e+12};

        std::vector<int> idx;
        for (auto ei : en) {
            idx.push_back(Numlib::round<int>(ei / egrain));
        }
        for (std::size_t i = 0; i < idx.size(); ++i) {
            CHECK(std::abs(wsum(idx[i]) - wsum_ans(i)) / wsum_ans(i) < 5.0e-2);
        }
    }

    SECTION("Ethane_hindered_rotor_steinrab") // Stein and Rabinovitch (1973)
    {
        double emax = 40000.0;
        double egrain = 5.0;
        int ngrains = 1 + Numlib::round<int>(emax / egrain);

        Numlib::Vec<double> vibr = {2915.0, 2915.0, 1388.0, 995.0,  1370.0,
                                    2974.0, 2974.0, 1460.0, 1460.0, 822.0,
                                    822.0,  2950.0, 2950.0, 1469.0, 1469.0,
                                    1190.0, 1190.0};

        double sigma = 3.0;
        double rotc = 10.704;
        double v0 = 1024.0;

        auto wsum = Sc::steinrab(vibr, sigma, rotc, v0, ngrains, egrain, true);

        Numlib::Vec<double> en = {500.0,  1000.0,  2000.0,  4000.0,  6000.0,
                                  8000.0, 10000.0, 20000.0, 30000.0, 40000.0};

        Numlib::Vec<double> wsum_ans = {
            2.0,       8.0,       55.7,      1187.0,     1.3111e+4,
            9.8719e+4, 5.7428e+5, 4.2130e+8, 4.4204e+10, 1.66923e+12};

        std::vector<int> idx;
        for (auto ei : en) {
            idx.push_back(Numlib::round<int>(ei / egrain));
        }
        for (std::size_t i = 0; i < idx.size(); ++i) {
            CHECK(std::abs(wsum(idx[i]) - wsum_ans(i)) / wsum_ans(i) < 0.17);
        }
    }
}

