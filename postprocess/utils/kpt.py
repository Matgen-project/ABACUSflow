#!/usr/bin/env python
# -*- coding: utf-8 -*-

from math import floor
import numpy as np

from .stru import Stru
import seekpath


class Kpt:
    def __init__(self, coords, lattice, weight=None, label=None, to_unit_cell=False,
                 coords_are_cartesian=False):
        self._lattice = lattice
        self._fcoords = self._lattice.frac_coords(coords) if coords_are_cartesian else coords
        self._label = label
        self._w = weight
        if to_unit_cell:
            for i in range(len(self._fcoords)):
                self._fcoords[i] -= floor(self._fcoords[i])
        self._ccoords = self._lattice.cart_coords(self._fcoords)

    @property
    def label(self):
        return self._label

    @property
    def frac_coords(self):
        return np.copy(self._fcoords)

    @property
    def cart_coords(self):
        return np.copy(self._ccoords)

    @property
    def a(self):
        return self._fcoords[0]

    @property
    def b(self):

        return self._fcoords[1]

    @property
    def c(self):
        return self._fcoords[2]

    @property
    def weight(self):
        return self._w

    def __str__(self):
        return "{} {}".format(self.frac_coords, self.label)

    @staticmethod
    def from_kpt(filename):
        with open(filename, "r") as f:
            lines = [i.strip() for i in f.readlines()]
        title, n, line_mode, *high_kpt = lines
        if '#' in line_mode:
            line_mode = line_mode.split('#')[0].strip(' ')
        if line_mode.lower() in ["direct", "cartesian"]:
            tmp = np.asarray(high_kpt.split(), dtype=int)
            density, shift = tmp[:3], tmp[3:]
            return density, shift, line_mode
        elif line_mode.lower in ["mp", "gamma"]:
            infos = []
            for i in high_kpt:
                if not i:
                    continue
                coord, lb = i.split('#')
                infos.append(
                    np.asarray(coord.split(), dtype=float)
                )
            tmp = np.asarray(infos)
            coords, weight = tmp[:, :3], tmp[:, -1]
            return coords, weight, line_mode, n
        elif line_mode.lower() == "line":
            infos, lbs = [], []
            for i in high_kpt:
                if not i:
                    continue
                coord, lb = i.split('#')
                infos.append(
                    np.asarray(coord.split(), dtype=float)
                )
                lbs.append(lb.strip(' '))
            tmp = np.asarray(infos)
            coords, number = tmp[:, :3], tmp[:, -1]

            return coords, number, lbs, line_mode, n
        else:
            raise RuntimeError("something wrong! check the KPT file!")

    @staticmethod
    def generate_kpath_from(filename, time_reversal=True):
        structure = Stru.from_stru(filename)
        cell = (
            structure.lattice.lattice,
            structure.frac_coords,
            structure.atomic_number
        )

        kpath = seekpath.get_explicit_k_path(cell, with_time_reversal=time_reversal)
        return kpath


if __name__ == "__main__":
    pass

