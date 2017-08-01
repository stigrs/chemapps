#!/usr/bin/env python

"""
Configure script for setting C/C++/Fortran compiler options in Makefiles.

@note Tested on Linux, Cygwin, Windows, and Mac OS X.

@author  Stig Rune Sellevag <stigrs@gmail.com>.
@warning No warranty offered.
@date    2012-2015
"""

import os
import os.path
import platform
import argparse

class Error(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class CC_compiler:
    def __init__(self, uname):
        self.uname            = uname
        self.list_linux       = ["gcc", "icc", "pathcc", "opencc"]
        self.list_cygwin      = ["gcc"]
        self.list_windows     = ["cl", "icl"]
        self.list_mac         = ["clang", "gcc"]
        self.default          = self._set_default()
        self.selected         = []
        self.debug            = []
        self.profile          = []
        self.cflags           = []
        self.cc_depflag       = []
        self.cc_shared_flag   = []
        self.cc_shared_ldflag = []
        self.cc_libexe        = []
        self.cc_linkexe       = []
        self.gcc    = {"depflag": "-MM", 
                       "langflag": "-ansi -pedantic",
                       "warnflag": "-Wall -Wextra -Wshadow",
                       "debugflag": "-g", 
                       "profileflag": "-O2 -pg", 
                       "optflag": "-O3",
                       "shared_flag": "-fPIC", 
                       "shared_ldflag": "-shared"}
        self.icc    = {"depflag": "-MM", 
                       "langflag": "-strict-ansi",
                       "warnflag": "-Wall",
                       "debugflag": "-g", 
                       "profileflag": "-O2 -pg", 
                       "optflag": "-O3 -ipo -no-prec-div -xHOST",
                       "shared_flag": "-fPIC", 
                       "shared_ldflag": "-shared"}     
        self.pathcc = {"depflag": "-MM", 
                       "langflag": "-ansi",
                       "warnflag": "-Wall -Wno-uninitialized -Wshadow",
                       "debugflag": "-fullwarn -g", 
                       "profileflag": "-O2 -pg", 
                       "optflag": "-Ofast",
                       "shared_flag": "-fPIC", 
                       "shared_ldflag": "-shared"}     
        self.opencc = {"depflag": "-MM", 
                       "langflag": "-ansi",
                       "warnflag": "-Wall -Wno-uninitialized -Wshadow",
                       "debugflag": "-fullwarn -g", 
                       "profileflag": "-O2 -pg", 
                       "optflag": "-Ofast",
                       "shared_flag": "-fPIC", 
                       "shared_ldflag": "-shared"}     
        self.clang  = {"depflag": "-MM", 
                       "langflag": "-pedantic",
                       "warnflag": "-Wall",
                       "debugflag": "-g", 
                       "profileflag": "-O2 -pg", 
                       "optflag": "-Ofast",
                       "shared_flag": "-fPIC", 
                       "shared_ldflag": "-shared"}
        self.cl  = {"depflag": "", 
                    "langflag": "/nologo /EHsc /Za",
                    "warnflag": 
                    "/D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS",
                    "debugflag": "/Zi", 
                    "profileflag": "", 
                    "optflag": "/O2",
                    "shared_flag": "", 
                    "shared_ldflag": "",
                    "libexe": "lib"}
        self.icl = {"depflag": "", 
                    "langflag": "/nologo /EHsc /Za",
                    "warnflag": 
                    "/W3 /D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS",
                    "debugflag": "/Zi /debug:all", 
                    "profileflag": "/Qprof-gen", 
                    "optflag": "/fast",
                    "shared_flag": "", 
                    "shared_ldflag": "",
                    "libexe": "xilib"}
        self.clang = {"depflag": "-MM", 
                      "langflag": "-ansi -pedantic",
                      "warnflag": "-Wall -Wextra -Wshadow",
                      "debugflag": "-g", 
                      "profileflag": "-O2 -pg", 
                      "optflag": "-O3",
                      "shared_flag": "-fPIC", 
                      "shared_ldflag": "-shared"}
 

    def _set_default(self):
        if self.uname == "Linux":
            return self.list_linux[0]
        elif self.uname == "CYGWIN":
            return self.list_cygwin[0]
        elif self.uname == "Windows":
            return self.list_windows[0]
        elif self.uname == "Darwin":
            return self.list_mac[0]
        else:
            return "gcc"

    def _set_flags(self):
        if self.selected == "gcc":
            self.cflags.append(self.gcc["langflag"])
            self.cflags.append(" ")
            self.cflags.append(self.gcc["warnflag"])
            self.cflags.append(" ")
            if self.debug:
                self.cflags.append(self.gcc["debugflag"])
                self.cflags.append(" ")
            elif self.profile:
                self.cflags.append(self.gcc["profileflag"])
                self.cflags.append(" ")
            else:
                self.cflags.append(self.gcc["optflag"])
                self.cflags.append(" ")
            self.cc_depflag       = self.gcc["depflag"]
            self.cc_shared_flag   = self.gcc["shared_flag"]
            self.cc_shared_ldflag = self.gcc["shared_ldflag"]
        elif self.selected == "icc":
            self.cflags.append(self.icc["langflag"])
            self.cflags.append(" ")
            self.cflags.append(self.icc["warnflag"])
            self.cflags.append(" ")
            if self.debug:
                self.cflags.append(self.icc["debugflag"])
                self.cflags.append(" ")
            elif self.profile:
                self.cflags.append(self.icc["profileflag"])
                self.cflags.append(" ")
            else:
                self.cflags.append(self.icc["optflag"])
                self.cflags.append(" ")
            self.cc_depflag       = self.icc["depflag"]
            self.cc_shared_flag   = self.icc["shared_flag"]
            self.cc_shared_ldflag = self.icc["shared_ldflag"]
        elif self.selected == "pathcc":
            self.cflags.append(self.pathcc["langflag"])
            self.cflags.append(" ")
            self.cflags.append(self.pathcc["warnflag"])
            self.cflags.append(" ")
            if self.debug:
                self.cflags.append(self.pathcc["debugflag"])
                self.cflags.append(" ")
            elif self.profile:
                self.cflags.append(self.pathcc["profileflag"])
                self.cflags.append(" ")
            else:
                self.cflags.append(self.pathcc["optflag"])
                self.cflags.append(" ")
            self.cc_depflag       = self.pathcc["depflag"]
            self.cc_shared_flag   = self.pathcc["shared_flag"]
            self.cc_shared_ldflag = self.pathcc["shared_ldflag"]
        elif self.selected == "opencc":
            self.cflags.append(self.opencc["langflag"])
            self.cflags.append(" ")
            self.cflags.append(self.opencc["warnflag"])
            self.cflags.append(" ")
            if self.debug:
                self.cflags.append(self.opencc["debugflag"])
                self.cflags.append(" ")
            elif self.profile:
                self.cflags.append(self.opencc["profileflag"])
                self.cflags.append(" ")
            else:
                self.cflags.append(self.opencc["optflag"])
                self.cflags.append(" ")
            self.cc_depflag       = self.opencc["depflag"]
            self.cc_shared_flag   = self.opencc["shared_flag"]
            self.cc_shared_ldflag = self.opencc["shared_ldflag"]
        elif self.selected == "clang":
            self.cflags.append(self.clang["langflag"])
            self.cflags.append(" ")
            self.cflags.append(self.clang["warnflag"])
            self.cflags.append(" ")
            if self.debug:
                self.cflags.append(self.clang["debugflag"])
                self.cflags.append(" ")
            elif self.profile:
                self.cflags.append(self.clang["profileflag"])
                self.cflags.append(" ")
            else:
                self.cflags.append(self.clang["optflag"])
                self.cflags.append(" ")
            self.cc_depflag       = self.clang["depflag"]
            self.cc_shared_flag   = self.clang["shared_flag"]
            self.cc_shared_ldflag = self.clang["shared_ldflag"]
        elif self.selected == "cl":
            self.cflags.append(self.cl["langflag"])
            self.cflags.append(" ")
            self.cflags.append(self.cl["warnflag"])
            self.cflags.append(" ")
            if self.debug:
                self.cflags.append(self.cl["debugflag"])
                self.cflags.append(" ")
            elif self.profile:
                self.cflags.append(self.cl["profileflag"])
                self.cflags.append(" ")
            else:
                self.cflags.append(self.cl["optflag"])
                self.cflags.append(" ")
            self.cc_libexe = self.cl["libexe"]
        elif self.selected == "icl":
            self.cflags.append(self.icl["langflag"])
            self.cflags.append(" ")
            self.cflags.append(self.icl["warnflag"])
            self.cflags.append(" ")
            if self.debug:
                self.cflags.append(self.icl["debugflag"])
                self.cflags.append(" ")
            elif self.profile:
                self.cflags.append(self.icl["profileflag"])
                self.cflags.append(" ")
            else:
                self.cflags.append(self.icl["optflag"])
                self.cflags.append(" ")
            self.cc_libexe = self.icl["libexe"]
        return "".join(self.cflags)

    def supported(self):
        if self.uname == "Linux":
            return self.list_linux
        elif self.uname == "CYGWIN":
            return self.list_cygwin
        elif self.uname == "Windows":
            return self.list_windows
        elif self.uname == "Darwin":
            return self.list_mac
        else:
            return self.default

    def get_default(self):
        return self.default

    def set(self, selected, debug, profile):
        self.selected = selected
        self.debug    = debug
        self.profile  = profile
        self.cflags   = self._set_flags()

    def name(self):
        return self.selected

    def flags(self):
        return self.cflags

    def depflag(self):
        return self.cc_depflag

    def shared_flag(self):
        return self.cc_shared_flag

    def shared_ldflag(self):
        return self.cc_shared_ldflag

    def libexe(self):
        return self.cc_libexe

class CXX_compiler:
    def __init__(self, uname):
        self.uname             = uname
        self.list_linux        = ["g++", "icpc", "pathCC", "openCC"]
        self.list_cygwin       = ["g++"]
        self.list_windows      = ["cl", "icl"]
        self.list_mac          = ["clang++", "g++"]
        self.default           = self._set_default()
        self.selected          = []
        self.debug             = []
        self.profile           = []
        self.cxxflags          = []
        self.cxx_depflag       = []
        self.cxx_shared_flag   = []
        self.cxx_shared_ldflag = []
        self.cxx_libexe        = []
        self.cxx_linkexe       = []
        self.gxx    = {"depflag": "-MM", 
                       "langflag": "-pedantic -std=c++11",
                       "warnflag": "-Wall -Wextra -Wshadow",
                       "debugflag": "-g", 
                       "profileflag": "-O2 -pg", 
                       "optflag": "-O3",
                       "shared_flag": "-fPIC", 
                       "shared_ldflag": "-shared"}
        self.icpc   = {"depflag": "-MM", 
                       "langflag": "-strict-ansi",
                       "warnflag": "-Wall",
                       "debugflag": "-g", 
                       "profileflag": "-O2 -pg", 
                       "optflag": "-O3 -ipo -no-prec-div -xHOST",
                       "shared_flag": "-fPIC", 
                       "shared_ldflag": "-shared"}     
        self.pathCC = {"depflag": "-MM", 
                       "langflag": "-std=c++98",
                       "warnflag": "-Wall -Wno-uninitialized -Wshadow",
                       "debugflag": "-fullwarn -g", 
                       "profileflag": "-O2 -pg", 
                       "optflag": "-Ofast",
                       "shared_flag": "-fPIC", 
                       "shared_ldflag": "-shared"}     
        self.openCC = {"depflag": "-MM", 
                       "langflag": "-std=c++98 -pedantic-errors",
                       "warnflag": "-Wall -Wextra -Wno-uninitialized -Wshadow",
                       "debugflag": "-g", 
                       "profileflag": "-O2 -pg", 
                       "optflag": "-Ofast",
                       "shared_flag": "-fPIC", 
                       "shared_ldflag": "-shared"}     
        self.clangxx = {"depflag": "-MM", 
                       "langflag": "-pedantic -std=c++11",
                       "warnflag": "-Wall",
                       "debugflag": "-g", 
                       "profileflag": "-O2 -pg", 
                       "optflag": "-Ofast",
                       "shared_flag": "-fPIC", 
                       "shared_ldflag": "-shared"}
        self.cl  = {"depflag": "", 
                    "langflag": "/nologo /EHsc /Za",
                    "warnflag": 
                    "/D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS",
                    "debugflag": "/Zi", 
                    "profileflag": "", 
                    "optflag": "/O2",
                    "shared_flag": "", 
                    "shared_ldflag": "",
                    "libexe": "lib"}
        self.icl = {"depflag": "", 
                    "langflag": "/nologo /EHsc /Za",
                    "warnflag": 
                    "/W3 /D_SCL_SECURE_NO_WARNINGS /D_CRT_SECURE_NO_WARNINGS",
                    "debugflag": "/Zi /debug:all", 
                    "profileflag": "/Qprof-gen", 
                    "optflag": "/fast",
                    "shared_flag": "", 
                    "shared_ldflag": "",
                    "libexe": "xilib"}

    def _set_default(self):
        if self.uname == "Linux":
            return self.list_linux[0]
        elif self.uname == "CYGWIN":
            return self.list_cygwin[0]
        elif self.uname == "Windows":
            return self.list_windows[0]
        elif self.uname == "Darwin":
            return self.list_mac[0]
        else:
            return "g++"

    def _set_flags(self):
        if self.selected == "g++":
            self.cxxflags.append(self.gxx["langflag"])
            self.cxxflags.append(" ")
            self.cxxflags.append(self.gxx["warnflag"])
            self.cxxflags.append(" ")
            if self.debug:
                self.cxxflags.append(self.gxx["debugflag"])
                self.cxxflags.append(" ")
            elif self.profile:
                self.cxxflags.append(self.gxx["profileflag"])
                self.cxxflags.append(" ")
            else:
                self.cxxflags.append(self.gxx["optflag"])
                self.cxxflags.append(" ")
            self.cxx_depflag       = self.gxx["depflag"]
            self.cxx_shared_flag   = self.gxx["shared_flag"]
            self.cxx_shared_ldflag = self.gxx["shared_ldflag"]
        elif self.selected == "icpc":
            self.cxxflags.append(self.icpc["langflag"])
            self.cxxflags.append(" ")
            self.cxxflags.append(self.icpc["warnflag"])
            self.cxxflags.append(" ")
            if self.debug:
                self.cxxflags.append(self.icpc["debugflag"])
                self.cxxflags.append(" ")
            elif self.profile:
                self.cxxflags.append(self.icpc["profileflag"])
                self.cxxflags.append(" ")
            else:
                self.cxxflags.append(self.icpc["optflag"])
                self.cxxflags.append(" ")
            self.cxx_depflag       = self.icpc["depflag"]
            self.cxx_shared_flag   = self.icpc["shared_flag"]
            self.cxx_shared_ldflag = self.icpc["shared_ldflag"]
        elif self.selected == "pathCC":
            self.cxxflags.append(self.pathCC["langflag"])
            self.cxxflags.append(" ")
            self.cxxflags.append(self.pathCC["warnflag"])
            self.cxxflags.append(" ")
            if self.debug:
                self.cxxflags.append(self.pathCC["debugflag"])
                self.cxxflags.append(" ")
            elif self.profile:
                self.cxxflags.append(self.pathCC["profileflag"])
                self.cxxflags.append(" ")
            else:
                self.cxxflags.append(self.pathCC["optflag"])
                self.cxxflags.append(" ")
            self.cxx_depflag       = self.pathCC["depflag"]
            self.cxx_shared_flag   = self.pathCC["shared_flag"]
            self.cxx_shared_ldflag = self.pathCC["shared_ldflag"]
        elif self.selected == "openCC":
            self.cxxflags.append(self.openCC["langflag"])
            self.cxxflags.append(" ")
            self.cxxflags.append(self.openCC["warnflag"])
            self.cxxflags.append(" ")
            if self.debug:
                self.cxxflags.append(self.openCC["debugflag"])
                self.cxxflags.append(" ")
            elif self.profile:
                self.cxxflags.append(self.openCC["profileflag"])
                self.cxxflags.append(" ")
            else:
                self.cxxflags.append(self.openCC["optflag"])
                self.cxxflags.append(" ")
            self.cxx_depflag       = self.openCC["depflag"]
            self.cxx_shared_flag   = self.openCC["shared_flag"]
            self.cxx_shared_ldflag = self.openCC["shared_ldflag"]
        elif self.selected == "clang++":
            self.cxxflags.append(self.clangxx["langflag"])
            self.cxxflags.append(" ")
            self.cxxflags.append(self.clangxx["warnflag"])
            self.cxxflags.append(" ")
            if self.debug:
                self.cxxflags.append(self.clangxx["debugflag"])
                self.cxxflags.append(" ")
            elif self.profile:
                self.cxxflags.append(self.clangxx["profileflag"])
                self.cxxflags.append(" ")
            else:
                self.cxxflags.append(self.clangxx["optflag"])
                self.cxxflags.append(" ")
            self.cxx_depflag       = self.clangxx["depflag"]
            self.cxx_shared_flag   = self.clangxx["shared_flag"]
            self.cxx_shared_ldflag = self.clangxx["shared_ldflag"]
        elif self.selected == "cl":
            self.cxxflags.append(self.cl["langflag"])
            self.cxxflags.append(" ")
            self.cxxflags.append(self.cl["warnflag"])
            self.cxxflags.append(" ")
            if self.debug:
                self.cxxflags.append(self.cl["debugflag"])
                self.cxxflags.append(" ")
            elif self.profile:
                self.cxxflags.append(self.cl["profileflag"])
                self.cxxflags.append(" ")
            else:
                self.cxxflags.append(self.cl["optflag"])
                self.cxxflags.append(" ")
            self.cxx_libexe = self.cl["libexe"]
        elif self.selected == "icl":
            self.cxxflags.append(self.icl["langflag"])
            self.cxxflags.append(" ")
            self.cxxflags.append(self.icl["warnflag"])
            self.cxxflags.append(" ")
            if self.debug:
                self.cxxflags.append(self.icl["debugflag"])
                self.cxxflags.append(" ")
            elif self.profile:
                self.cxxflags.append(self.icl["profileflag"])
                self.cxxflags.append(" ")
            else:
                self.cxxflags.append(self.icl["optflag"])
                self.cxxflags.append(" ")
            self.cxx_libexe = self.icl["libexe"]
        return "".join(self.cxxflags)

    def supported(self):
        if self.uname == "Linux":
            return self.list_linux
        elif self.uname == "CYGWIN":
            return self.list_cygwin
        elif self.uname == "Windows":
            return self.list_windows
        elif self.uname == "Darwin":
            return self.list_mac
        else:
            return self.default

    def get_default(self):
        return self.default

    def set(self, selected, debug, profile):
        self.selected = selected
        self.debug    = debug
        self.profile  = profile
        self.cxxflags = self._set_flags()

    def name(self):
        return self.selected

    def flags(self):
        return self.cxxflags

    def depflag(self):
        return self.cxx_depflag

    def shared_flag(self):
        return self.cxx_shared_flag

    def shared_ldflag(self):
        return self.cxx_shared_ldflag

    def libexe(self):
        return self.cxx_libexe

class FC_compiler:
    def __init__(self, uname):
        self.uname            = uname
        self.list_linux       = ["gfortran", "ifort", "pathf95", "openf95"]
        self.list_cygwin      = ["gfortran"]
        self.list_windows     = ["ifort"]
        self.list_mac         = ["gfortran"]
        self.default          = self._set_default()
        self.selected         = []
        self.debug            = []
        self.profile          = []
        self.fcflags          = []
        self.fc_depflag       = []
        self.fc_shared_flag   = []
        self.fc_shared_ldflag = []
        self.gfortran = {"modflag": "-J", 
                         "langflag": "-std=f95 -fall-intrinsics -fdefault-integer-8",
                         "warnflag": "-Wall -Wextra",
                         "debugflag": "-g", 
                         "profileflag": "-pg", 
                         "optflag": "-O3",
                         "shared_flag": "-fPIC", 
                         "shared_ldflag": "-shared"}
        self.ifort    = {"modflag": "-module ", 
                         "langflag": "-std95 -diag-disable 7416 -i8",
                         "warnflag": "-warn",
                         "debugflag": "-g", 
                         "profileflag": "-pg", 
                         "optflag": "-O3 -ipo -no-prec-div -xHOST",
                         "shared_flag": "-fPIC", 
                         "shared_ldflag": "-shared"}   
        self.ifortw   = {"modflag": "/module:", 
                         "langflag": "/stand:f95 /Qdiag-disable:7416",
                         "warnflag": "/warn:all",
                         "debugflag": "/Zi /debug:all", 
                         "profileflag": "/Qprof-gen", 
                         "optflag": "/fast",
                         "shared_flag": "", 
                         "shared_ldflag": ""}           
        self.pathf95  = {"modflag": "-module ", 
                         "langflag": "",
                         "warnflag": "-Wall",
                         "debugflag": "-g", 
                         "profileflag": "-pg", 
                         "optflag": "-Ofast",
                         "shared_flag": "-fPIC", 
                         "shared_ldflag": "-shared"}
        self.openf95  = {"modflag": "-module ", 
                         "langflag": "",
                         "warnflag": "-Wall",
                         "debugflag": "-g", 
                         "profileflag": "-pg", 
                         "optflag": "-Ofast",
                         "shared_flag": "-fPIC", 
                         "shared_ldflag": "-shared"}

    def _set_default(self):
        if self.uname == "Linux":
            return self.list_linux[0]
        elif self.uname == "CYGWIN":
            return self.list_cygwin[0]
        elif self.uname == "Windows":
            return self.list_windows[0]
        elif self.uname == "Darwin":
            return self.list_mac[0]
        else:
            return "gfortran"

    def _set_flags(self):
        if self.selected == "gfortran":
            self.fcflags.append(self.gfortran["langflag"])
            self.fcflags.append(" ")
            self.fcflags.append(self.gfortran["warnflag"])
            self.fcflags.append(" ")
            if self.debug:
                self.fcflags.append(self.gfortran["debugflag"])
                self.fcflags.append(" ")
            elif self.profile:
                self.fcflags.append(self.gfortran["profileflag"])
                self.fcflags.append(" ")
            else:
                self.fcflags.append(self.gfortran["optflag"])
                self.fcflags.append(" ")
            self.fc_modflag       = self.gfortran["modflag"]
            self.fc_shared_flag   = self.gfortran["shared_flag"]
            self.fc_shared_ldflag = self.gfortran["shared_ldflag"]
        elif self.selected == "ifort":
            self.fcflags.append(self.ifort["langflag"])
            self.fcflags.append(" ")
            self.fcflags.append(self.ifort["warnflag"])
            self.fcflags.append(" ")
            if self.debug:
                self.fcflags.append(self.ifort["debugflag"])
                self.fcflags.append(" ")
            elif self.profile:
                self.fcflags.append(self.ifort["profileflag"])
                self.fcflags.append(" ")
            else:
                self.fcflags.append(self.ifort["optflag"])
                self.fcflags.append(" ")
            self.fc_modflag       = self.ifort["modflag"]
            self.fc_shared_flag   = self.ifort["shared_flag"]
            self.fc_shared_ldflag = self.ifort["shared_ldflag"]
        elif self.selected == "ifortw":
            self.fcflags.append(self.ifortw["langflag"])
            self.fcflags.append(" ")
            self.fcflags.append(self.ifortw["warnflag"])
            self.fcflags.append(" ")
            if self.debug:
                self.fcflags.append(self.ifortw["debugflag"])
                self.fcflags.append(" ")
            elif self.profile:
                self.fcflags.append(self.ifortw["profileflag"])
                self.fcflags.append(" ")
            else:
                self.fcflags.append(self.ifortw["optflag"])
                self.fcflags.append(" ")
            self.fc_modflag       = self.ifortw["modflag"]
            self.fc_shared_flag   = self.ifortw["shared_flag"]
            self.fc_shared_ldflag = self.ifortw["shared_ldflag"]
        elif self.selected == "pathf95":
            self.fcflags.append(self.pathf95["langflag"])
            self.fcflags.append(" ")
            self.fcflags.append(self.pathf95["warnflag"])
            self.fcflags.append(" ")
            if self.debug:
                self.fcflags.append(self.pathf95["debugflag"])
                self.fcflags.append(" ")
            elif self.profile:
                self.fcflags.append(self.pathf95["profileflag"])
                self.fcflags.append(" ")
            else:
                self.fcflags.append(self.pathf95["optflag"])
                self.fcflags.append(" ")
            self.fc_modflag       = self.pathf95["modflag"]
            self.fc_shared_flag   = self.pathf95["shared_flag"]
            self.fc_shared_ldflag = self.pathf95["shared_ldflag"]
        elif self.selected == "openf95":
            self.fcflags.append(self.openf95["langflag"])
            self.fcflags.append(" ")
            self.fcflags.append(self.openf95["warnflag"])
            self.fcflags.append(" ")
            if self.debug:
                self.fcflags.append(self.openf95["debugflag"])
                self.fcflags.append(" ")
            elif self.profile:
                self.fcflags.append(self.openf95["profileflag"])
                self.fcflags.append(" ")
            else:
                self.fcflags.append(self.openf95["optflag"])
                self.fcflags.append(" ")
            self.fc_modflag       = self.openf95["modflag"]
            self.fc_shared_flag   = self.openf95["shared_flag"]
            self.fc_shared_ldflag = self.openf95["shared_ldflag"]
        return "".join(self.fcflags)

    def supported(self):
        if self.uname == "Linux":
            return self.list_linux
        elif self.uname == "CYGWIN":
            return self.list_cygwin
        elif self.uname == "Windows":
            return self.list_windows
        elif self.uname == "Darwin":
            return self.list_mac
        else:
            return self.default

    def get_default(self):
        return self.default

    def set(self, selected, debug, profile):
        self.selected = selected
        if (self.uname == "Windows") and (selected == "ifort"):
            self.selected = "ifortw"
        self.debug   = debug
        self.profile = profile
        self.fcflags = self._set_flags()

    def name(self):
        if self.selected == "ifortw":
            return "ifort"
        else:
            return self.selected

    def flags(self):
        return self.fcflags

    def modflag(self):
        return self.fc_modflag

    def shared_flag(self):
        return self.fc_shared_flag

    def shared_ldflag(self):
        return self.fc_shared_ldflag

def get_platform():
    """Get operating system."""
    uname = platform.system()
    if uname.find("CYGWIN") != -1:
        uname = "CYGWIN"
    return uname

def get_exeext(uname):
    """Get extension of executables."""
    if (uname == "Windows") or (uname == "CYGWIN"):
        return ".exe"
    else:
        return ""

def which(program):
    """Function that mimics the behavior of the UNIX 'which' command."""
    def _is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    def _ext_candidates(fpath):
        yield fpath
        for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
            yield fpath + ext

    fpath, fname = os.path.split(program)
    if fpath:
        if _is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            for candidate in _ext_candidates(exe_file):
                if _is_exe(candidate):
                    return candidate
    return None

def list_supported(cc, cxx, fc):
    """List supported compilers."""
    print("Supported compiler(s): ")
    print("C:       ", cc.supported())
    print("C++:     ", cxx.supported())
    print("Fortran: ", fc.supported())

def find_compiler(comp_selected, comp_list):
    """Find selected compiler in list of supported compilers and in PATH."""
    found = False
    for c in comp_list:
        if c == comp_selected:
            found = True
    if not found:
        raise Error(comp_selected)
    if not which(comp_selected):
        raise Error(comp_selected)        

def find_makefiles(uname):
    template = "Makefile.unix"
    if uname == "Windows":
        template = "Makefile.win"
    currdir  = os.getcwd()
    filelist = []
    for root, dirs, files in os.walk(os.getcwd()):
        for name in files:
            if name == template:
                filelist.append(os.path.relpath(os.path.join(root, name)))
    return filelist

def generate_makefile(filename, uname, debug, profile, prefix, libdir, 
                      includedir, cc, cxx, fc):
    exeext   = get_exeext(uname)

    defs = ""
    if uname == "Windows":
        ldflags  = "/LIBDIR:" + "\"" + libdir + "\""
        cppflags = "/I" + "\"" + includedir + "\""
        if debug:
            defs = ""
        else:
            defs = "/DARMA_NO_DEBUG"
    else:
        ldflags  = "-L" + "\"" + libdir + "\""
        cppflags = "-I" + "\"" + includedir + "\""
        if debug:
            defs = ""
        else:
            defs = "-DARMA_NO_DEBUG"


    data = []
    finp = open(filename, "r")
    for line in finp:
        if line.find("@prefix@") != -1:
            data.append(line.replace("@prefix@", prefix))
        elif line.find("@cc@") != -1:
            data.append(line.replace("@cc@", cc.name()))
        elif line.find("@ccflags@") != -1:
            data.append(line.replace("@ccflags@", cc.flags()))
        elif line.find("@cc_depflag@") != -1:
            data.append(line.replace("@cc_depflag@", cc.depflag()))
        elif line.find("@cc_shared_flag@") != -1:
            data.append(line.replace("@cc_shared_flag@", cc.shared_flag()))
        elif line.find("@cc_shared_ldflag@") != -1:
            data.append(line.replace("@cc_shared_ldflag@", cc.shared_ldflag()))
        elif line.find("@cc_libexe@") != -1:
            data.append(line.replace("@cc_libexe@", cc.libexe()))
        elif line.find("@cxx@") != -1:
            data.append(line.replace("@cxx@", cxx.name()))
        elif line.find("@cxxflags@") != -1:
            data.append(line.replace("@cxxflags@", cxx.flags()))
        elif line.find("@cxx_depflag@") != -1:
            data.append(line.replace("@cxx_depflag@", cxx.depflag()))
        elif line.find("@cxx_shared_flag@") != -1:
            data.append(line.replace("@cxx_shared_flag@", cxx.shared_flag()))
        elif line.find("@cxx_shared_ldflag@") != -1:
            data.append(line.replace("@cxx_shared_ldflag@",cxx.shared_ldflag()))
        elif line.find("@cxx_libexe@") != -1:
            data.append(line.replace("@cxx_libexe@", cxx.libexe()))
        elif line.find("@fc@") != -1:
            data.append(line.replace("@fc@", fc.name()))
        elif line.find("@fcflags@") != -1:
            data.append(line.replace("@fcflags@", fc.flags()))
        elif line.find("@fc_depflag@") != -1:
            data.append(line.replace("@fc_depflag@", fc.depflag()))
        elif line.find("@fc_modflag@") != -1:
            data.append(line.replace("@fc_modflag@", fc.modflag()))
        elif line.find("@fc_shared_flag@") != -1:
            data.append(line.replace("@fc_shared_flag@", fc.shared_flag()))
        elif line.find("@fc_shared_ldflag@") != -1:
            data.append(line.replace("@fc_shared_ldflag@", fc.shared_ldflag()))
        elif line.find("@defs@") != -1:
            data.append(line.replace("@defs@", defs))
        elif line.find("@cppflags@") != -1:
            data.append(line.replace("@cppflags@", cppflags))
        elif line.find("@ldflags@") != -1:
            data.append(line.replace("@ldflags@", ldflags))
        elif line.find("@exeext@") != -1:
            data.append(line.replace("@exeext@", exeext))
        else:
            data.append(line)
    finp.close()

    head, tile = os.path.split(filename)
    makefile   = os.path.join(head, "Makefile")

    print("Configuring ", makefile, " ...")

    fout = open(makefile, "w")
    for line in data:
        fout.write(line)
    fout.close()
            
def configure():
    """This is the driver of the configure script."""

    # Determine operating system and set some default values:

    uname = get_platform()
    home  = "HOME"
    if uname == "Windows":
        home = "USERPROFILE"

    default_prefix     = os.environ[home]
    default_libdir     = os.path.join(default_prefix, "lib")
    default_includedir = os.path.join(default_prefix, "include")

    # Initialize compilers:

    cc  = CC_compiler(uname)
    cxx = CXX_compiler(uname)
    fc  = FC_compiler(uname)

    # Parse input arguments:

    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--list", action="store_true", dest="list",
                        default=False, help="list supported compilers and exit")
    parser.add_argument("-d", "--debug", action="store_true", dest="debug", 
                        default=False, help="perform debugging")
    parser.add_argument("-p", "--profile", action="store_true", dest="profile",
                        default=False, help="perform profiling")
    parser.add_argument("--prefix", action="store", dest="prefix", 
                        default=default_prefix, 
                        help="installation directory [$HOME/bin]")
    parser.add_argument("--libdir", action="store", dest="libdir", 
                        default=default_libdir, 
                        help="installation directory [$HOME/lib]")
    parser.add_argument("--includedir", action="store", dest="includedir", 
                        default=default_includedir, 
                        help="installation directory [$HOME/include]")
    parser.add_argument("--CC", action="store", dest="CC", 
                        default=cc.get_default(), help="C compiler")
    parser.add_argument("--CXX", action="store", dest="CXX", 
                        default=cxx.get_default(), help="C++ compiler")
    parser.add_argument("--FC", action="store", dest="FC", 
                        default=fc.get_default(), help="Fortran compiler")
    args = parser.parse_args()

    prefix     = args.prefix
    libdir     = args.libdir
    includedir = args.includedir
    debug      = args.debug
    profile    = args.profile

    if debug and profile:
        debug = False

    if args.list:
        list_supported(cc, cxx, fc)
        return

    # Check if selected compiler is supported and set compiler flags:

    try:
        find_compiler(args.CC, cc.supported())
        cc.set(args.CC, debug, profile)
    except Error as e:
        raise Error("Unsupported C compiler: " + e.value)
    try:
        find_compiler(args.CXX, cxx.supported())
        cxx.set(args.CXX, debug, profile)
    except Error as e:
        raise Error("Unsupported C++ compiler: " + e.value)
    try:
        find_compiler(args.FC, fc.supported())
        fc.set(args.FC, debug, profile)
    except Error as e:
        raise Error("Unsupported Fortran compiler: " + e.value)

    # Output information:

    print("Operating system:    ", uname)
    print("")
    print("Installation prefix: ", prefix)
    print("")
    print("C compiler:          ", cc.name(), cc.flags())
    print("C++ compiler:        ", cxx.name(), cxx.flags())
    print("Fortran compiler:    ", fc.name(), fc.flags())
    print("")
    print("libdir:              ", libdir)
    print("includedir:          ", includedir)
    print("")

    if debug:
        print("Debugging is specified\n")
    if profile:
        print("Profiling is specified\n")

    # Create Makefiles:

    makefiles = find_makefiles(uname)
    if len(makefiles) == 0:
        raise Error("no makefile templates found")
    for f in makefiles:
        generate_makefile(f, uname, debug, profile, prefix, 
                          libdir, includedir, cc, cxx, fc)

    make = "make"
    if uname == "Windows":
        make = "nmake"

    print("\nYou now have the following options:\n")
    print("  %s           - compile program/library"      % make)
    print("  %s tests     - compile/run test cases"       % make)
    print("  %s doc       - compile documentation"        % make)
    print("  %s install   - install program/library"      % make)
    print("  %s uninstall - unistall program/library"     % make)
    print("  %s clean     - remove some generated files"  % make)
    print("  %s distclean - remove all generated files\n" % make)
    
if __name__ == "__main__":
    try:
        configure()
    except Error as e:
        print(e.value)

