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
#include <srs/datum.h>
#include <srs/math.h>
#include <boost/math/special_functions/gamma.hpp>
#include <cmath>


srs::dvector statecount::bswine(const srs::dvector& vibr,
                                int ngrains,
                                double egrain,
                                bool sum,
                                const srs::dvector& rot)
{
    srs::dvector result = srs::zeros<srs::dvector>(ngrains);
    if (!rot.empty()) {  // initialize with rotational states
        result = rot;
    }
    else {
        if (sum) {  // count sum of states
            result = srs::ones<srs::dvector>(ngrains);
        }
        else {  // count density of states
            result(0) = 1.0;
        }
    }
    for (auto w : vibr) {
        int wj = srs::round<int>(w / egrain);
        for (int i = wj; i < ngrains; ++i) {
            result(i) += result(i - wj);
        }
    }
    if (!sum) {
        result *= 1.0 / egrain;
    }
    return result;
}
srs::dvector statecount::free_rotor(
    int sigma, double rotc, int ngrains, double egrain, bool sum)
{
    srs::dvector result = srs::zeros<srs::dvector>(ngrains);

    int k    = 0;
    double f = -0.5;

    if (sum) {
        k = 1;
        f = 0.5;
    }
    double qr = std::sqrt(datum::pi) / (sigma * std::sqrt(rotc));
    double gf = boost::math::tgamma<double>(f + 1.0);

    for (int i = 0; i < ngrains; ++i) {
        result(i) = qr * std::pow(i * egrain, f) / gf;
    }
    return result;
}
