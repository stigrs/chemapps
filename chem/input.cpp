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

#include <chem/arma_io.h>
#include <chem/input.h>
#include <chem/ptable.h>
#include <chem/utils.h>
#include <gsl/gsl>

std::istream& operator>>(std::istream& from, Input& inp)
{
    switch (inp.type) {
    case Input::t_int:
        inp.read_int(from);
        break;
    case Input::t_long:
        inp.read_long(from);
        break;
    case Input::t_uint:
        inp.read_uint(from);
        break;
    case Input::t_ulint:
        inp.read_ulint(from);
        break;
    case Input::t_double:
        inp.read_double(from);
        break;
    case Input::t_string:
        inp.read_string(from);
        break;
    case Input::t_ivector:
        inp.read_ivector(from);
        break;
    case Input::t_uvector:
        inp.read_uvector(from);
        break;
    case Input::t_dvector:
        inp.read_dvector(from);
        break;
    case Input::t_mol_formula:
        inp.read_mol_formula(from);
        break;
    case Input::t_noval:
        throw Input_invalid("cannot read data with no type");
    default:
        throw Input_invalid("cannot read data with unknown type");
    }
    inp.state = Input::init;
    return from;
}

std::ostream& operator<<(std::ostream& to, const Input& inp)
{
    if (inp.state == Input::not_init) {
        to << "not initialized";
    }
    else {
        arma::ivec iv;
        arma::uvec uv;
        arma::vec dv;
        std::vector<Mol_formula> mf;

        chem::Format<char> line;
        chem::Format<double> fix;
        line.width(32).fill('-');
        fix.fixed().width(8).precision(4);

        switch (inp.type) {
        case Input::t_int:
            to << *static_cast<int*>(inp.data);
            break;
        case Input::t_long:
            to << *static_cast<long*>(inp.data);
            break;
        case Input::t_uint:
            to << *static_cast<unsigned*>(inp.data);
            break;
        case Input::t_ulint:
            to << *static_cast<unsigned long*>(inp.data);
            break;
        case Input::t_double:
            to << *static_cast<double*>(inp.data);
            break;
        case Input::t_string:
            to << *static_cast<std::string*>(inp.data);
            break;
        case Input::t_ivector:
            iv = *static_cast<arma::ivec*>(inp.data);
            chem::print_vector(to, iv);
            break;
        case Input::t_uvector:
            uv = *static_cast<arma::uvec*>(inp.data);
            chem::print_vector(to, uv);
            break;
        case Input::t_dvector:
            dv = *static_cast<arma::vec*>(inp.data);
            chem::print_vector(to, dv);
            break;
        case Input::t_mol_formula:
            mf = *static_cast<std::vector<Mol_formula>*>(inp.data);
            to << line('-') << '\n'
               << "Center\tAtomic\tStoich.\tMass/amu\n"
               << "Number\tSymbol\n"
               << line('-') << '\n';
            for (std::size_t i = 0; i < mf.size(); ++i) {
                to << i + 1 << '\t' << mf[i].atom << '\t' << mf[i].stoich
                   << '\t' << fix(ptable::get_atomic_mass(mf[i].atom)) << '\n';
            }
            to << line('-') << '\n';
            break;
        default:
            throw Input_invalid("unknown type; cannot write data");
        }
    }
    return to;
}

//-----------------------------------------------------------------------------

void Input::read_int(std::istream& from)
{
    int& i = *static_cast<int*>(data);

    std::string buf;
    from >> buf;

    i = chem::from_string<int>(buf);
}

void Input::read_long(std::istream& from)
{
    long& i = *static_cast<long*>(data);

    std::string buf;
    from >> buf;

    i = chem::from_string<long>(buf);
}

void Input::read_uint(std::istream& from)
{
    unsigned& i = *static_cast<unsigned*>(data);

    std::string buf;
    from >> buf;

    i = chem::from_string<unsigned>(buf);
}

void Input::read_ulint(std::istream& from)
{
    unsigned long& i = *static_cast<unsigned long*>(data);

    std::string buf;
    from >> buf;

    i = chem::from_string<unsigned long>(buf);
}

void Input::read_double(std::istream& from)
{
    double& d = *static_cast<double*>(data);

    std::string buf;
    from >> buf;

    d = chem::from_string<double>(buf);
}

void Input::read_string(std::istream& from)
{
    std::string& s = *static_cast<std::string*>(data);

    std::string buf;
    from >> buf;

    s = chem::from_string<std::string>(buf);
}

void Input::read_ivector(std::istream& from)
{
    arma::ivec& v = *static_cast<arma::ivec*>(data);
    chem::read_vector(from, v);
}

void Input::read_uvector(std::istream& from)
{
    arma::uvec& v = *static_cast<arma::uvec*>(data);
    chem::read_vector(from, v);
}

void Input::read_dvector(std::istream& from)
{
    arma::vec& v = *static_cast<arma::vec*>(data);
    chem::read_vector(from, v);
}

void Input::read_mol_formula(std::istream& from)
{
    std::vector<Mol_formula>& v = *static_cast<std::vector<Mol_formula>*>(data);

    std::string buf;
    from >> buf;

    int n = chem::from_string<int>(buf);
    Expects(n >= 1);
    v.resize(n);

    char ch;
    std::string atom;
    int stoich;

    from >> ch;
    if (ch != '[') {
        throw Input_IO_error("'[' missing in molecular formula");
    }
    for (int i = 0; i < n; ++i) {
        from >> atom >> stoich >> ch;
        if (!from) {
            throw Input_IO_error("found no data for molecular formula");
        }
        if ((ch != ',') && (ch != ';') && (ch != ']')) {
            throw Input_IO_error("bad separator in molecular formula: "
                                 + chem::to_string(ch));
        }
        if (ch == ']') {
            from.unget();
        }
        if (!ptable::atomic_symbol_is_valid(atom)) {
            throw Input_invalid("bad atomic symbol: " + atom);
        }
        if (stoich < 1) {
            throw Input_invalid("bad stoichiometry: "
                                + chem::to_string(stoich));
        }
        v[i].atom   = atom;
        v[i].stoich = stoich;
    }
    from >> ch;
    if (ch != ']') {
        throw Input_IO_error("']' missing in molecular formula");
    }
}