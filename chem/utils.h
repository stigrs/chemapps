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

#ifndef CHEM_UTILS_H
#define CHEM_UTILS_H

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>

namespace chem {

// Assertion methods:

// Stroustrup's templated Assert() function (p. 751 in TC++PL).
template <typename E, typename A>
inline void Assert(A assertion)
{
    if (!assertion) {
        throw E();
    }
}

// Stroustrup's templated Assert() function (p. 752 in TC++PL).
template <typename A, typename E>
inline void Assert(A assertion, E except)
{
    if (!assertion) {
        throw except;
    }
}

// Stream handling methods:

// Derived class for reporting file opening errors.
struct Fopen_error : std::runtime_error {
    Fopen_error(std::string s) : std::runtime_error(s) {}
};

// Derived class for reporting lexical cast errors.
struct Bad_lexical_cast : std::bad_cast {
    const char* what() const throw() { return "bad lexical cast"; }
};

// Wrapper function for opening file stream for input.
inline void fopen(std::ifstream& from,
                  const std::string& filename,
                  std::ios_base::openmode mode = std::ios_base::in)
{
    from.open(filename.c_str(), mode);
    if (!from.is_open()) {
        throw Fopen_error("cannot open " + filename);
    }
}

// Wrapper function for opening file stream for output.
inline void fopen(std::ofstream& to,
                  const std::string& filename,
                  std::ios_base::openmode mode = std::ios_base::out)
{
    to.open(filename.c_str(), mode);
    if (!to.is_open()) {
        throw Fopen_error("cannot open " + filename);
    }
}

// Find section in input stream.
inline bool find_section(std::istream& from, const std::string& key)
{
    from.clear();
    from.seekg(0, std::ios_base::beg);

    std::string buf;
    while (std::getline(from, buf)) {
        std::istringstream iss(buf);
        iss >> buf;
        if (buf == key) {
            return true;
        }
    }
    return false;
}

// Lexical cast (type must be able to stream into and/or out of a string).
template <typename Target, typename Source>
Target lexical_cast(Source arg)
{
    std::stringstream interpreter;
    Target result;

    if (!(interpreter << arg) || !(interpreter >> result) ||
        !(interpreter >> std::ws).eof()) {  // stuff left in stream?
        throw Bad_lexical_cast();
    }
    return result;
}

// Stream format methods:

template <typename T>
struct Bound_form;

template <typename T>
std::ostream& operator<<(std::ostream& os, const Bound_form<T>& bf);

// Stroustrup's format class, pp. 635-636 in TC++PL (slightly modified).
template <typename T>
class Format {
public:
    explicit Format(int p = 6)
        : prc(p), wdt(0), ch(' '), fmt(std::ios_base::fmtflags(0))
    {
    }

    Bound_form<T> operator()(T vv) const { return Bound_form<T>(*this, vv); }

    // Set width.
    Format<T>& width(int w)
    {
        wdt = w;
        return *this;
    }

    // Set fill character.
    Format<T>& fill(char c)
    {
        ch = c;
        return *this;
    }

    // Set precision.
    Format<T>& precision(int p)
    {
        prc = p;
        return *this;
    }

    // Set general format.
    Format<T>& general()
    {
        fmt = std::ios_base::fmtflags(0);
        return *this;
    }

    // Set fixed floating point format.
    Format<T>& fixed()
    {
        fmt = std::ios_base::fixed;
        return *this;
    }

    // Set scientific floating point format using e character.
    Format<T>& scientific()
    {
        fmt = std::ios_base::scientific;
        return *this;
    }

    // Set scientific floating point format using E character.
    Format<T>& scientific_E()
    {
        fmt = std::ios_base::scientific | std::ios_base::uppercase;
        return *this;
    }

private:
    int prc;  // precision, default precision is 6
    int wdt;  // width, 0 means as wide as necessary
    char ch;  // fill character, default is blank

    std::ios_base::fmtflags fmt;  // format flag value

    friend std::ostream& operator<<<>(std::ostream&, const Bound_form<T>&);
};  // Format

// Form plus value.
template <typename T>
struct Bound_form {
    const Format<T>& f;

    T v;

    Bound_form(const Format<T>& ff, T vv) : f(ff), v(vv) {}

private:
    Bound_form& operator=(const Bound_form&);
};

template <typename T>
std::ostream& operator<<(std::ostream& to, const Bound_form<T>& bf)
{
    std::ostringstream s;

    s.width(bf.f.wdt);
    s.precision(bf.f.prc);
    s.fill(bf.f.ch);

    s << std::setiosflags(bf.f.fmt) << bf.v;
    return to << s.str();
}

// String methods:

// Derived class for reporting string cast errors.
struct Bad_from_string : std::runtime_error {
    Bad_from_string(std::string s) : std::runtime_error(s) {}
};

// Derived class for handling string find errors.
struct String_find_error : std::runtime_error {
    String_find_error(std::string s) : runtime_error(s) {}
};

// Simple convert to string method.
template <typename T>
std::string to_string(const T& t)
{
    std::ostringstream oss;
    oss << t;
    return oss.str();
}

// Simple extract from string method.
template <typename T>
T from_string(const std::string& s)
{
    std::istringstream iss(s);
    T t;
    if (!(iss >> t)) {
        throw Bad_from_string("bad cast from string '" + s + "'");
    }
    return t;
}

// Convert Fortran scientific D format to double.
inline double from_fortran_sci_fmt(const std::string& s)
{
    std::string ss             = s;
    std::string::size_type pos = ss.find('D');
    if (pos != std::string::npos) {
        ss.replace(pos, 1, 1, 'e');
    }
    return from_string<double>(ss);
}

// Trim leading and trailing characters from string.
inline std::string trim(const std::string& str, const char* sep)
{
    const std::string::size_type pos = str.find_first_not_of(sep);
    return (pos == std::string::npos)
               ? std::string()
               : str.substr(pos, str.find_last_not_of(sep) - pos + 1);
}

// Strip suffix from filename.
inline std::string strip_suffix(const std::string& filename,
                                const std::string& suffix)
{
    std::string basename       = filename;
    std::string::size_type pos = basename.rfind(suffix);
    if (pos == std::string::npos) {
        throw String_find_error(filename + " does not contain " + suffix);
    }
    return basename.erase(pos, basename.size() - pos);
}

// Get suffix from filename.
inline std::string get_suffix(const std::string& filename)
{
    std::string suffix         = filename;
    std::string::size_type pos = suffix.rfind(".");
    if (pos == std::string::npos) {
        throw String_find_error(filename + " does not have a suffix");
    }
    return suffix.erase(0, pos);
}

// Check if string has only blank characters.
inline bool str_has_only_blanks(const std::string& str)
{
    if (!str.empty()) {
        if (str.find_first_not_of(" \t") != std::string::npos) {
            return false;
        }
    }
    return true;
}

// Case-insensitive string comparison.
//
// Note: Can only handle characters in the default C locale!
inline bool stricmp(const std::string& str1, const std::string& str2)
{
    std::string s1;
    std::string s2;
    std::string::const_iterator it;
    for (it = str1.begin(); it < str1.end(); ++it) {
        if (::isalpha(*it)) {  // only convert alphabetic characters
            s1 += static_cast<char>(::tolower(*it));
        }
        else {
            s1 += *it;
        }
    }
    for (it = str2.begin(); it < str2.end(); ++it) {
        if (::isalpha(*it)) {
            s2 += static_cast<char>(::tolower(*it));
        }
        else {
            s2 += *it;
        }
    }
    return s1 == s2;
}
}  // namespace chem

#endif  // CHEM_UTILS_H
