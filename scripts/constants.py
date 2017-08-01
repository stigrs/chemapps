#   This file is part of PyChem - Python Chemistry Toolkit
#
#   Copyright (C) 2016-2017  Stig Rune Sellevag
#
#   PyChem is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   PyChem is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Module providing physical constants and metric prefixes."""

import os
import sys
from pychem.lib.chem.quantity import Quantity

# _CONSTANTS is defined outside of class as a look-up table to avoid
# time-consuming reading of the file with element data every time
# an instance is created.

_CONSTANTS = {}


def get(label):
    "Get physical constant."""
    return _CONSTANTS[label].value


def get_wunit(label):
    """Get physical constant with unit."""
    return _CONSTANTS[label]


def list_constants(fname=sys.stdout):
    """List all physical constants to output file object."""
    for key, value in sorted(_CONSTANTS.items()):
        fname.write("%-55s = %s\n" % (key, value))


def _load_constants():
    """Load physical constants from data file."""
    filename = os.path.join("pychem", "data", "constants.txt")
    with open(filename, "r") as finp:
        for line in finp:
            if line[0] == "#":
                continue
            else:
                data = line.split("\t")
                _CONSTANTS[data[0]] = Quantity(
                    float(data[1]), data[3][:-1], fmt="%.9e")


def _init():
    """Initialize table of physical constants."""
    if not _CONSTANTS:
        _load_constants()  # this is slow

# Initialize table of physical constants when the module is imported:
_init()

# Metric prefixes:

yotta = 1.0e+24
zetta = 1.0e+21
exa = 1.0e+18
peta = 1.0e+15
tera = 1.0e+12
giga = 1.0e+9
mega = 1.0e+6
kilo = 1.0e+3
hecto = 1.0e+2
deca = 10.0
one = 1.0
deci = 1.0e-1
centi = 1.0e-2
milli = 1.0e-3
micro = 1.0e-6
nano = 1.0e-9
pico = 1.0e-12
femto = 1.0e-15
atto = 1.0e-18
zepto = 1.0e-21
yocto = 1.0e-24

# Commonly used physical constants:

mu = get("atomic mass constant")
NA = get("Avogadro constant")
a0 = get("Bohr radius") * 1.0e+10  # 1E-10 m
k = get("Boltzmann constant")
G0 = get("conductance quantum")
epsilon0 = get("electric constant")
me = get("electron mass")
eV = get("electron volt")
e = get("elementary charge")
F = get("Faraday constant")
alpha = get("fine-structure constant")
R = get("molar gas constant")
Eh = get("Hartree energy")
c = get("speed of light in vacuum")
mu0 = get("mag. constant")
Phi0 = get("mag. flux quantum")
mp_me = get("proton-electron mass ratio")
G = get("Newtonian constant of gravitation")
h = get("Planck constant")
hbar = get("Planck constant over 2 pi")
mp = get("proton mass")
Rinf = get("Rydberg constant")
atm = get("standard atmosphere")
sigma = get("Stefan-Boltzmann constant")
