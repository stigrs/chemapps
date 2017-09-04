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
#include <armadillo>
#include <catch/catch.hpp>
#include <iostream>

TEST_CASE("test_math")
{
    SECTION("even_odd")
    {
        int even_number = 4;
        CHECK(chem::is_even(even_number) == true);

        int odd_number = -3;
        CHECK(chem::is_odd(odd_number) == true);
    }

    SECTION("integration")
    {
        // Abramowitz and Stegun, Table 25.4, p. 916:
        arma::vec xans = {-0.960289856497536,
                          -0.796666477413627,
                          -0.525532409916329,
                          -0.183434642495659,
                          0.183434642495659,
                          0.525532409916329,
                          0.796666477413627,
                          0.960289856497536};

        // Abramowitz and Stegun, Table 25.4, p. 916:
        arma::vec wans = {0.101228536290376,
                          0.222381034453374,
                          0.313706645877887,
                          0.362683783378362,
                          0.362683783378362,
                          0.313706645877887,
                          0.222381034453374,
                          0.101228536290376};

        int n = 8;
        arma::vec x(n);
        arma::vec w(n);
        chem::gaussleg(n, x, w);

        CHECK(arma::approx_equal(x, xans, "absdiff", 1.0e-10));
        CHECK(arma::approx_equal(w, wans, "absdiff", 1.0e-10));

        // Numpy:
        xans = {1.01985507,
                1.10166676,
                1.2372338,
                1.40828268,
                1.59171732,
                1.7627662,
                1.89833324,
                1.98014493};

        // Numpy:
        wans = {0.05061427,
                0.11119052,
                0.15685332,
                0.18134189,
                0.18134189,
                0.15685332,
                0.11119052,
                0.05061427};

        chem::gaussleg(n, x, w, 1.0, 2.0);

        CHECK(arma::approx_equal(x, xans, "absdiff", 1.0e-8));
        CHECK(arma::approx_equal(w, wans, "absdiff", 1.0e-8));
    }
}
