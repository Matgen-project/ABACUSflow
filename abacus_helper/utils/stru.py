#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np

from ._base import Latt, order_table


class Stru:
    def __init__(self, formula_symbol, formula_cell, elements, formula_positions, atom_number, species_number=None,
                 atom_potential=None, atom_basis=None, is_cart=False):
        self._fs = formula_symbol
        self._e = elements
        if species_number is None:
            species_number = []
            for item in self._e:
                species_number.append(
                    int(order_table.get(item))
                )
        self._sn = species_number
        self._a = atom_number
        self._lattice = Latt(formula_cell)
        self._fcoords = self._lattice.frac_coords(formula_positions) if is_cart else formula_positions
        self._ap = atom_potential
        self._ab = atom_basis
        self._ccoords = self._lattice.cart_coords(self._fcoords)

    @property
    def formula(self):
        return self._fs

    @property
    def lattice(self):
        return self._lattice

    @property
    def frac_coords(self):
        return self._fcoords

    @property
    def cart_coords(self):
        return self._ccoords

    @property
    def atoms_base(self):
        return self._ab

    @property
    def atom_potential(self):
        return self._ap

    @property
    def atom_species(self):
        return self._e

    @property
    def uniq_atom_species(self):
        uniq = []
        for i in self._e:
            if i not in uniq:
                uniq.append(i)

        return uniq

    @property
    def atom_number(self):
        return self._a

    @property
    def atomic_number(self):
        return self._sn

    @classmethod
    def from_stru(cls, filename='STRU'):
        # Read structure information from abacus structure file
        with open(filename, "r") as stru:
            lines = stru.readlines()
        # initialize reading information
        temp = []
        for line in lines:
            line = line.strip()
            line = line.replace('\n', ' ')
            line = line.replace('\t', ' ')
            line = line.replace('//', ' ')
            line = line.replace('#', ' ')

            if len(line) != 0:
                temp.append(line)

        atom_species = 0
        for i in range(len(temp)):
            if temp[i] == 'NUMERICAL_ORBITAL':
                atom_species = i - 1
                break

        atom_symbol = []
        atom_mass = []
        atom_potential = []
        atom_number = []
        atom_magnetism = []
        atom_positions = []
        atom_fix = []

        # get symbol, mass, potential
        for i in range(1, atom_species + 1):
            atom_symbol.append(temp[i].split()[0])
            atom_mass.append(float(temp[i].split()[1]))
            atom_potential.append(temp[i].split()[2])
            atom_number.append(0)
            atom_magnetism.append(0)
            atom_positions.append([])
            atom_fix.append([])

        # get basis
        atom_basis = []
        for i in range(atom_species + 2, (atom_species + 1) * 2):
            atom_basis.append(temp[i].split()[0])

        # get lattice
        atom_lattice_scale = float(temp[(atom_species + 1) * 2 + 1].split()[0])
        atom_lattice = np.array(
            [[float(temp[(atom_species + 1) * 2 + 3 + i].split()[:3][j])
              for j in range(3)] for i in range(3)])

        # get coordinates type
        atom_coor = temp[(atom_species + 1) * 2 + 7].split()[0]

        # get position,  atoms number, magnetism, fix
        for i in range(atom_species):
            pos_start = (atom_species + 1) * 2 + 8 + 3 * i
            for j in range(i):
                pos_start += atom_number[j]
            atom_it = atom_symbol.index(temp[pos_start].split()[0])
            atom_magnetism[atom_it] = float(temp[pos_start + 1].split()[0])
            atom_number[atom_it] = int(temp[pos_start + 2].split()[0])

            atom_positions[atom_it] = np.array(
                [[float(temp[pos_start + 3 + i].split()[:3][j])
                  for j in range(3)] for i in range(atom_number[atom_it])])

            atom_fix[atom_it] = np.array(
                [[int(temp[pos_start + 3 + i].split()[3:6][j])
                  for j in range(3)] for i in range(atom_number[atom_it])])

        # Reset structure information and return results
        formula_symbol = ''
        formula_positions = []
        for i in range(atom_species):
            if atom_number[i] == 1:
                formula_symbol += atom_symbol[i]

            else:
                formula_symbol += atom_symbol[i] + str(atom_number[i])

            for j in range(atom_number[i]):
                formula_positions.append(atom_positions[i][j])

        formula_cell = atom_lattice * atom_lattice_scale * 0.529177210903

        if atom_coor == 'Direct':
            is_cart = False

        elif atom_coor == 'Cartesian':
            is_cart = True
        else:
            raise ValueError("atomic coordinate type is ERROR")
        species = []
        for i, e in enumerate(atom_symbol):
            species.extend(
                [e, ] * atom_number[i]
            )

        return cls(formula_symbol, formula_cell, species, np.asarray(formula_positions),
                   atom_number, atom_potential=atom_potential,
                   atom_basis=atom_basis, is_cart=is_cart)


if __name__ == "__main__":
    pass


