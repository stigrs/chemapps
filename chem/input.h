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

#ifndef CHEM_INPUT_H
#define CHEM_INPUT_H

#include <armadillo>
#include <iostream>
#include <stdexcept>
#include <string>

// Error reporting:

struct Input_invalid : std::invalid_argument {
    Input_invalid(std::string s) : std::invalid_argument(s) {}
};

struct Input_IO_error : std::runtime_error {
    Input_IO_error(std::string s) : std::runtime_error(s) {}
};

// Forward declarations to allow friend declarations:

class Input;

std::istream& operator>>(std::istream& from, Input& inp);
std::ostream& operator<<(std::ostream& to, const Input& inp);

//
// Class for reading input data into a map.
//
class Input {
public:
    Input() : data(0), type(t_noval), state(not_init) {}

    explicit Input(int& i) : data(&i), type(t_int), state(not_init) {}
    explicit Input(long& i) : data(&i), type(t_long), state(not_init) {}
    explicit Input(unsigned& i) : data(&i), type(t_uint), state(not_init) {}
    explicit Input(unsigned long& i) : data(&i), type(t_ulint), state(not_init)
    {
    }

    Input(int& i, int ii) : data(&i), type(t_int), state(def) { i = ii; }
    Input(long& i, long ii) : data(&i), type(t_long), state(def) { i = ii; }
    Input(unsigned& i, unsigned ii) : data(&i), type(t_uint), state(def)
    {
        i = ii;
    }

    Input(unsigned long& i, unsigned long ii)
        : data(&i), type(t_ulint), state(def)
    {
        i = ii;
    }

    Input(double& d) : data(&d), type(t_double), state(not_init) {}
    Input(double& d, double dd) : data(&d), type(t_double), state(def)
    {
        d = dd;
    }

    Input(std::string& s) : data(&s), type(t_string), state(not_init) {}
    Input(std::string& s, std::string ss) : data(&s), type(t_string), state(def)
    {
        s = ss;
    }

    Input(arma::ivec& v) : data(&v), type(t_ivector), state(not_init) {}
    Input(arma::ivec& v, arma::ivec& vv) : data(&v), type(t_ivector), state(def)
    {
        v = vv;
    }

    Input(arma::uvec& v) : data(&v), type(t_uvector), state(not_init) {}
    Input(arma::uvec& v, arma::uvec& vv) : data(&v), type(t_uvector), state(def)
    {
        v = vv;
    }

    Input(arma::vec& v) : data(&v), type(t_dvector), state(not_init) {}
    Input(arma::vec& v, arma::vec& vv) : data(&v), type(t_dvector), state(def)
    {
        v = vv;
    }

    bool is_init() const { return !(state == not_init); }
    bool is_default() const { return state == def; }

    friend std::istream& operator>>(std::istream& from, Input& inp);
    friend std::ostream& operator<<(std::ostream& to, const Input& inp);

private:
    /// Definition of valid types.
    enum Type {
        t_noval,    // no value
        t_int,      // int
        t_long,     // long
        t_uint,     // unsigned int
        t_ulint,    // unsigned long int
        t_double,   // double
        t_string,   // std::string
        t_ivector,  // integer vector
        t_uvector,  // unsigned vector
        t_dvector   // double vector
    };

    /// Definition of intialization states.
    enum State { init, not_init, def };

    void* data;
    Type type;
    State state;

    void read_int(std::istream& from);
    void read_long(std::istream& from);
    void read_uint(std::istream& from);
    void read_ulint(std::istream& from);
    void read_double(std::istream& from);
    void read_string(std::istream& from);
    void read_ivector(std::istream& from);
    void read_uvector(std::istream& from);
    void read_dvector(std::istream& from);
};  // Input

#endif  // CHEM_INPUT_H
