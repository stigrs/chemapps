/**
   @file argparse.cpp
   
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

#include <chem/argparse.h>
#include <chem/utils.h>
#include <cstring>

//----------------------------------------------------------------------------

int    Argparse::argc;
char** Argparse::argv;
char   Argparse::invalid[20];

//----------------------------------------------------------------------------

void Argparse::init(int argc_, char** argv_)
{
    std::strcpy(invalid, "____invalid____");
    argc = argc_;
    argv = argv_;
}

bool Argparse::has_switch(const char* arg)
{
    const char* this_arg;

    int    tmp_argc = argc;
    char** tmp_argv = argv;

    while (++tmp_argv && --tmp_argc) {
	this_arg = *tmp_argv;
	if (! std::strcmp(this_arg, arg) && (tmp_argc > 0)) {
	    return true;
	}
    }
    return false;
}

const char* Argparse::read(const char* arg, const char* def)
{
    const char* this_arg;
    
    int    tmp_argc = argc;
    char** tmp_argv = argv;

    while (++tmp_argv && --tmp_argc) {
	this_arg = *tmp_argv;
	if (! std::strcmp(this_arg, arg) && (tmp_argc > 0)) {
	    return *++tmp_argv;
	}
    }
    return def; // no match, return default value
}

double Argparse::read(const char* arg, const double& def)
{
    const char* pc = read(arg, invalid);

    if (std::strcmp(pc, invalid) == 0) {
	return def;
    }
    else {
	return chem::from_string<double>(pc);
    }
}

int Argparse::read(const char* arg, const int& def)
{
    const char* pc = read(arg, invalid);

    if (std::strcmp(pc, invalid) == 0) {
	return def;
    }
    else {
	return chem::from_string<int>(pc);
    }
}

std::string Argparse::read(const char* arg, const std::string& def)
{
    const char* pc = read(arg, invalid);

    if (std::strcmp(pc, invalid) == 0) {
	return def;
    }
    else {
	return pc;
    }
}
