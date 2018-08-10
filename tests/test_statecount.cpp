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

#include <chem/statecount.h>
#include <srs/array.h>
#include <srs/math.h>
#include <catch/catch.hpp>
#include <iostream>

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
}
