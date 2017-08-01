/**
   @file argparse.h
   
   This file is part of ChemApps - A C++ Chemistry Toolkit
   
   Copyright (C) 2016-2017  Stig Rune Sellevag
   
   ChemApps is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   ChemApps is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CHEM_ARGPARSE_H
#define CHEM_ARGPARSE_H

#include <iostream>
#include <string>


/**
   Class for reading command line arguments with or without value.
   
   Based on code written by Roger Hansen <rogerha@ifi.uio.no> and
   Hans Petter Langtangen <hpl@ifi.uio.no>.
   
   @note Static only.
*/
class Argparse {
public:
    static void init(int argc_, char** argv_);

    static bool has_switch(const char* arg);

    static std::string read(const char* arg, const std::string& def);
    static const char* read(const char* arg, const char* def);
    static double      read(const char* arg, const double& def);
    static int         read(const char* arg, const int& def);

private:
    static int    argc;
    static char** argv;
    static char   invalid[20];
};

#endif /* CHEM_ARGPARSE_H */
