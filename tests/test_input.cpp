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

#include <chem/input.h>
#include <chem/utils.h>
#include <armadillo>
#include <catch/catch.hpp>
#include <map>
#include <string>

TEST_CASE("test_input")
{
    int i;
    double d;
    std::string s;

    arma::ivec iv;
    arma::uvec uv;
    arma::vec dv;
    arma::ivec iv_ans = {1, 2, 3, 4};
    arma::uvec uv_ans = {5, 6, 7, 8};
    arma::vec dv_ans  = {0.1, 0.2, 0.3, 0.4, 0.5};

    std::map<std::string, Input> data;
    data["integer"] = Input(i);
    data["double"]  = Input(d);
    data["string"]  = Input(s);
    data["ivector"] = Input(iv);
    data["uvector"] = Input(uv);
    data["dvector"] = Input(dv);

    std::ifstream from;
    chem::fopen(from, "test_input.inp");

    std::string key;
    while (from >> key) {
        auto it = data.find(key);
        if (it != data.end()) {
            from >> it->second;
        }
    }

    CHECK(i == 1);
    CHECK(d == 2.0);
    CHECK(s == "hello");

    for (arma::uword it = 0; it < iv_ans.size(); ++it) {
        CHECK(iv(it) == iv_ans(it));
    }
    for (arma::uword it = 0; it < uv_ans.size(); ++it) {
        CHECK(uv(it) == uv_ans(it));
    }
    for (arma::uword it = 0; it < dv_ans.size(); ++it) {
        CHECK(dv(it) == dv_ans(it));
    }
}
