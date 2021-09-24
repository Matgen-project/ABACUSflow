#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Modules for handing crystallographic lattice-parameters.
"""

import numpy as np


class Latt:
    """
    Construct Lattice parameter object.
    """

    def __init__(self, matrix=None):
        if not isinstance(matrix, np.ndarray):
            _m = np.asarray(matrix,
                            dtype=np.float64).reshape((3, 3))
        else:
            _m = matrix
        self._lat = _m.round(decimals=5)
        self._a, self._b, self._c = self._lat

    def __repr__(self):
        return "{}\n{}\n{}".format(*self.lattice)

    @property
    def a(self):
        return np.linalg.norm(self._a)

    @property
    def b(self):
        return np.linalg.norm(self._b)

    @property
    def c(self):
        return np.linalg.norm(self._c)

    @property
    def va(self):
        return self._a

    @property
    def vb(self):
        return self._b

    @property
    def vc(self):
        return self._c

    @classmethod
    def from_parameters(cls, a, b, c, alpha, beta, gamma):
        """
        Construct a new Lattice object from parameters.
        Args:
            a:
            b:
            c:
            alpha:
            beta:
            gamma:

        Returns:

        """
        angles_r = np.radians([alpha, beta, gamma])
        cos_alpha, cos_beta, cos_gamma = np.cos(angles_r)
        sin_alpha, sin_beta, sin_gamma = np.sin(angles_r)
        val = cls._abs_cap((cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta))
        va = [a * sin_beta, 0.0, a * cos_beta]
        vb = [-b * sin_alpha * np.cos(np.arccos(val)),
              b * sin_alpha * np.sin(np.arccos(val)), b * cos_alpha]
        vc = [0.0, 0.0, float(c)]
        return cls(np.asarray([va, vb, vc]))

    @staticmethod
    def _abs_cap(val, max_abs_val=1):
        """
        Return the value with its absolute value capped at max_abs_val.

        Particularly useful in passing values to trignometric functions where
        numerical errors may result in an argument > 1 being passed in.

        Args:

            val (float): Input value.

            max_abs_val (float): The maximum absolute value for val. Defaults to 1.

        Returns:
            val if abs(val) < 1 else sign of val * max_abs_val.
        """
        return max(min(val, max_abs_val), -max_abs_val)

    @classmethod
    def cubic(cls, a):
        """Construct cubic Lattice from lattice parameter information."""
        return cls.from_parameters(a, a, a, 90, 90, 90)

    @classmethod
    def tetragonal(cls, a, c):
        """Construct tetragonal Lattice from lattice parameter information."""
        return cls.from_parameters(a, a, c, 90, 90, 90)

    @classmethod
    def orthorhombic(cls, a, b, c):
        """Construct orthorhombic Lattice."""
        return cls.from_parameters(a, b, c, 90, 90, 90)

    @classmethod
    def monoclinic(cls, a, b, c, beta):
        """Construct monoclinic Lattice from lattice parameter information."""
        return cls.from_parameters(a, b, c, 90, beta, 90)

    @classmethod
    def hexagonal(cls, a, c):
        """Construct hexagonal Lattice from lattice parameter information."""
        return cls.from_parameters(a, a, c, 90, 90, 120)

    @classmethod
    def rhombohedral(cls, a, alpha):
        """Construct rhombohedral Lattice."""
        return cls.from_parameters(a, a, a, alpha, alpha, alpha)

    @property
    def lattice(self):
        """

        Returns:  lattice matrix.

        """
        return self._lat

    def inv_lattice(self):
        """

        Returns: inverse lattice matrix.

        """
        return np.linalg.inv(self._lat)

    def reciprocal_lattice(self):
        """Return reciprocal Lattice."""
        return Latt(2 * np.pi * np.linalg.inv(self._lat).T)

    def reciprocal_lattice_crystallographic(self):
        """Return reciprocal Lattice without 2 * pi."""
        return Latt(self.reciprocal_lattice().lattice / (2 * np.pi))

    def divide_kpts(self, threshold=0.02):
        if threshold == 0:
            return [1, ] * 3
        rlc = self.reciprocal_lattice_crystallographic()
        ratio = []
        for i in [rlc.a, rlc.b, rlc.c]:
            ratio.append(np.floor(i / threshold))
        return np.asarray(ratio).astype(int).tolist()

    def draw_kpath(self):
        pass

    def cart_coords(self, frac_coords):
        """

        Args:
            frac_coords: fraction coords

        Returns: cartesian coords from fractional coords using Lattice.

        """
        return np.dot(np.array(frac_coords), self._lat)

    def frac_coords(self, cart_coords):
        """

        Args:
            cart_coords: cartesian coords

        Returns: fractional coords from cartesian coords using Lattice.

        """
        return np.dot(np.array(cart_coords), self.inv_lattice())


if __name__ == "__main__":
    pass
