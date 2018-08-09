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
#include <srs/math.h>
#include <srs/utils.h>
#include <catch/catch.hpp>
#include <fstream>
#include <string>


TEST_CASE("test_gaussnmr")
{
    std::ifstream from;
    srs::fopen(from, "test_gaussnmr.inp");

    std::string nmr_method = "MP2 GIAO";
    double degen_tol       = 0.05;

    std::vector<Gauss_NMR> nmr;
    std::vector<double> shield;

    Gauss_data gauss(from, out);
    gauss.get_nmr_data(nmr, nmr_method, degen_tol);

    for (auto& i : nmr) {
        shield.push_back(std::accumulate(i.shield.begin(), i.shield.end(), 0.0)
                         / static_cast<double>(i.shield.size()));
    }
    CHECK(shield.size() == 3);
    CHECK(srs::approx_equal(shield[0], 31.6691, 1.0e-12));
    CHECK(srs::approx_equal(shield[1], 197.6891, 1.0e-12));
    CHECK(srs::approx_equal(shield[2], 362.2317, 1.0e-12));
}
