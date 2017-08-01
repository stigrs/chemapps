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

"""
Module providing the periodic table of elements.
"""

import os
import sys
import re
from element import Element

# _ISOTOPES is defined outside of class as a look-up table to avoid
# time-consuming reading of the file with element data every time
# an instance is created.
_ISOTOPES = []


def get(number=None, symbol=""):
    """
    Return element from the periodic table of elements.

    The most abundant isotope is returned for a given atomic number or
    symbol.
    """
    elements = []
    for elem in _ISOTOPES:
        if number == elem.atomic_number():
            elements.append(elem)
        elif symbol == elem.atomic_symbol():
            elements.append(elem)
        elif symbol == elem.isotope_symbol():
            elements.append(elem)

    atom = elements[0]
    if len(elements) > 1:
        most_abundant = atom.isotope_comp()
        for elem in elements:
            if elem.isotope_comp() > most_abundant:
                most_abundant = elem.isotope_comp()
                atom = elem
    return atom


def list_isotopes():
    """Return all isotopes in the periodic table."""
    elements = []
    for elem in _ISOTOPES:
        elements.append(elem.isotope_symbol())
    return elements


def list_elements():
    """Return all elements in the periodic table."""
    elements = []
    for elem in _ISOTOPES:
        elements.append(elem.atomic_symbol())
    return [ii for n, ii in enumerate(elements) if ii not in elements[:n]]


def _load_isotopes():
    """Load isotopes from data file."""
    filename = "elements.txt"
    table = []
    isotope = {}
    with open(filename, "r") as finp:
        for line in finp:
            if line[0] == "#":
                continue
            elif line and (not line.isspace()):
                key, val = line.split("=")
                isotope[key] = val
            else:
                table.append(isotope)
                isotope = {}

    element = {}
    pattern1 = r"((\d+(\.\d*)?|\.\d+)?\((.*)?\)?)"
    pattern2 = r"\[(\d+(\.\d*)?|\.\d+)?,(\d+(\.\d*)?|\.\d+)?\]"
    pattern3 = r"\[(\d+)?\]"
    for item in table:
        at_number = int(item["Atomic Number "].strip(" \n"))
        mass_number = int(item["Mass Number "].strip(" \n"))
        symbol = item["Atomic Symbol "].strip(" \n")
        if symbol == "D" or symbol == "T":
            aze = symbol
        else:
            aze = str(mass_number) + symbol
        mass = item["Relative Atomic Mass "].strip(" \n")
        match = re.search(pattern1, mass)
        if match:
            mass = float(match.group(2))
        weight = item["Standard Atomic Weight "].strip(" \n")
        match1 = re.search(pattern1, weight)
        match2 = re.search(pattern2, weight)
        match3 = re.search(pattern3, weight)
        if match1:
            weight = float(match1.group(2))
        elif match2:
            wa = float(match2.group(1))
            wb = float(match2.group(3))
            weight = (wa + wb) / 2.0
        elif match3:
            weight = float(match3.group(1))
        else:
            weight = 0.0
        isocomp = item["Isotopic Composition "].strip(" \n")
        match = re.search(pattern1, isocomp)
        if match:
            isocomp = float(match.group(2))
        else:
            isocomp = 0.0
        element["Atomic Number"] = at_number
        element["Atomic Symbol"] = symbol
        element["Isotope Symbol"] = aze
        element["Mass Number"] = mass_number
        element["Relative Atomic Mass"] = mass
        element["Standard Atomic Weight"] = weight
        element["Isotopic Composition"] = isocomp
        elem = Element(element)
        _ISOTOPES.append(elem)
        element = {}


def _init():
    """Initialize periodic table of elements."""
    if not _ISOTOPES:
        _load_isotopes()  # this is slow

# Initialize periodic table of elements when the module is imported:
_init()

if __name__ == "__main__":
    fout = sys.stdout

    isotopes = list_isotopes()
    atoms = list_elements()

    fout.write("namespace ptable {\n")
    fout.write("    const int num_isotopes = %d\n" % len(isotopes))
    fout.write("    const int num_atoms = %d\n" % len(atoms))

    it = 0
    fout.write("    const std::string isotopes[num_isotopes] = {\n")
    fout.write("        ")
    for iso in isotopes:
        fout.write("\"%s\", " % iso)
        if it == 7:
            fout.write("\n        ");
            it = 0
        it += 1
    fout.write(" };\n");

    it = 0
    fout.write("    const int mass_numbers[num_isotopes] = {\n")
    fout.write("        ")
    for iso in isotopes:
        fout.write("%d, " % get(symbol=iso).mass_number())
        if it == 10:
            fout.write("\n        ");
            it = 0
        it += 1
    fout.write(" };\n");

    it = 0
    fout.write("    const double atomic_masses[num_isotopes] = {\n")
    fout.write("        ")
    for iso in isotopes:
        fout.write("%f, " % get(symbol=iso).atomic_mass())
        if it == 5:
            fout.write("\n        ");
            it = 0
        it += 1
    fout.write(" };\n");

    it = 0
    fout.write("    const double atomic_weights[num_isotopes] = {\n")
    fout.write("        ")
    for iso in isotopes:
        fout.write("%f, " % get(symbol=iso).atomic_weight())
        if it == 5:
            fout.write("\n        ");
            it = 0
        it += 1
    fout.write(" };\n");

    it = 0
    fout.write("    const double isotope_comp[num_isotopes] = {\n")
    fout.write("        ")
    for iso in isotopes:
        fout.write("%f, " % get(symbol=iso).isotope_comp())
        if it == 5:
            fout.write("\n        ");
            it = 0
        it += 1
    fout.write(" };\n");

    it = 0
    fout.write("    const double isotope_comp[num_isotopes] = {\n")
    fout.write("        ")
    for atom in atoms:
        fout.write("%d-%s, " % (get(symbol=atom).mass_number(),
                                get(symbol=atom)))
        if it == 5:
            fout.write("\n        ");
            it = 0
        it += 1
    fout.write(" };\n");

    it = 0
    fout.write("    const int atomic_numbers[num_isotopes] = {\n")
    fout.write("        ")
    for iso in isotopes:
        fout.write("%d, " % (get(symbol=iso).atomic_number()))
        if it == 10:
            fout.write("\n        ");
            it = 0
        it += 1
    fout.write(" };\n");
