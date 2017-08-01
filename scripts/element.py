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

"""Module providing class for representation of elements."""


class Element(object):
    """
    Class providing representation of an element in the periodic table.

    Attributes:
        Atomic Number
        Atomic Symbol
        Isotope Symbol
        Mass Number
        Relative Atomic Mass
        Standard Atomic Weight
        Isotopic Composition
    """
    def __init__(self, data):
        self.__data = data

    def __str__(self):
        """Return a human-readable string representation of the object."""
        return self.isotope_symbol()

    def __repr__(self):
        """Return a representation that can be used to represent the object."""
        text = """Atomic Number = %s
Atomic Symbol = %s
Mass Number = %s
Relative Atomic Mass = %s
Standard Atomic Weight = %s
Isotopic Composition = %s
"""
        return text % (self.atomic_number(), self.atomic_symbol(),
                       self.mass_number(), self.atomic_mass(),
                       self.atomic_weight(), self.isotope_comp())

    def get_data(self):
        """Get all data."""
        return self.__data

    def atomic_number(self):
        """Return atomic number."""
        return self.__data["Atomic Number"]

    def atomic_symbol(self):
        """Return atomic symbol."""
        return self.__data["Atomic Symbol"]

    def isotope_symbol(self):
        """Return isotope symbol."""
        return self.__data["Isotope Symbol"]

    def mass_number(self):
        """Return mass number."""
        return self.__data["Mass Number"]

    def atomic_mass(self):
        """Return relative atomic mass."""
        return self.__data["Relative Atomic Mass"]

    def atomic_weight(self):
        """Return standard atomic weight."""
        return self.__data["Standard Atomic Weight"]

    def isotope_comp(self):
        """Return isotopic composition."""
        return self.__data["Isotopic Composition"]
