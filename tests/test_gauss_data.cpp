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

#include <chem/gauss_data.h>
#include <srs/array.h>
#include <srs/math.h>
#include <srs/utils.h>
#include <catch/catch.hpp>
#include <fstream>


TEST_CASE("test_gauss_data")
{
    std::ifstream from;
    srs::fopen(from, "test_gauss_data_h2o.inp");

    srs::dvector ans
        = {-1.89865925E-04, -1.76751570E-16, 8.04584647E-01,  -3.41586331E-16,
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
    srs::dvector hess;

    Gauss_data gauss(from, fchk);
    gauss.get_hessians(hess);

    CHECK(hess.size() == ans.size());
    CHECK(srs::approx_equal(hess, ans, 1.0e-12));
}
