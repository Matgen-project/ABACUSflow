#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import numpy as np

from monty.io import reverse_readfile

from ._base import Spin


def read_statinfo(filename):
    with open(filename, "r") as info:
        lines = [i.strip() for i in info.readlines()]
    rgx = re.compile(r"\s\s+")
    titles = rgx.split(lines[0])
    if len(titles) == 7:
        ispin = True
    else:
        ispin = False
    _energy, kpoints, occupation = [], [], []
    nkpt, _n = 0, 0
    for line in lines:
        if line.startswith("BAND"):
            head_lst = rgx.split(line)
            kpt = np.asarray(
                head_lst[-1].strip('(').strip(')').split(), dtype=float
            )
            kpoints.append(kpt)
            nkpt += 1
        else:
            if line:
                _n += 1
                lst = rgx.split(line)
                _energy.append(
                    np.asarray(lst, dtype=float)
                )

    nbands = int(_n / nkpt)
    energy = [_energy[i: i + nbands] for i in range(0, _n, nbands)]

    return _merge(energy, ispin, nkpt, nbands)


def _merge(energy, ispin, nkpt, nbands, up_index=(1, 2), down_index=(3, 4)):
    if ispin:
        eigenvalues = {
            Spin.up: np.zeros((nkpt, nbands, 2)),
            Spin.down: np.zeros((nkpt, nbands, 2)),
        }
    else:
        eigenvalues = {Spin.up: np.zeros((nkpt, nbands, 2))}

    ikpt = -1
    a, b = up_index
    c, d = down_index
    for sl in energy:
        ikpt += 1
        for i in range(nbands):
            tmp = sl[i]
            if len(tmp) == 3:
                eigenvalues[Spin.up][ikpt, i, 0] = tmp[a]
                eigenvalues[Spin.up][ikpt, i, 1] = tmp[b]
            else:
                eigenvalues[Spin.up][ikpt, i, 0] = tmp[a]
                eigenvalues[Spin.up][ikpt, i, 1] = tmp[b]
                eigenvalues[Spin.down][ikpt, i, 0] = tmp[c]
                eigenvalues[Spin.down][ikpt, i, 1] = tmp[d]

    return eigenvalues


def read_pdos(filename):
    energy_values, energy_flag = [], False
    orbital_values, orbital_flag = [], False
    atom_pdos, atom_pdos_flag = [], False
    index, atom_index, species, l, m, z = [None, ] * 6

    def clean(eqs):
        _, val = eqs.split('=')
        res = "".join(val.split()).strip("\"")
        return int(res) if res.isdigit() else res

    pdos_dat = []
    with open(filename, "r") as f:
        for i in f:
            if "<energy_values units=\"eV\">" in i:
                energy_flag = True
                continue
            if "</energy_values>" in i:
                energy_flag = False
            if energy_flag:
                energy_values.append(i.strip())
            if "<orbital" in i:
                orbital_flag = True
                continue
            if "</orbital>" in i:
                orbital_flag = False
            if orbital_flag:
                if "index" in i and "atom" not in i:
                    index = clean(i)
                if "atom_index" in i:
                    atom_index = clean(i)
                if "species" in i:
                    species = clean(i)
                if 'l' in i:
                    l = clean(i)
                if 'm' in i:
                    m = clean(i)
                if 'z' in i:
                    z = clean(i)
                if "<data>" in i:
                    atom_pdos_flag = True
                    continue
                if "</data>" in i:
                    atom_pdos_flag = False
                    ipdos = {
                        "index": index, "atom_index": atom_index, "species": species,
                        "l": l, "m": m, "z": z, "pdos": np.asarray(atom_pdos, dtype=float)
                    }
                    pdos_dat.append(ipdos)
                    atom_pdos = []
                if atom_pdos_flag:
                    atom_pdos.append(i.strip().split())
    return pdos_dat, np.asarray(energy_values, dtype=float)


def read_orbital(filename):
    dats = []
    with open(filename, "r") as orbital:
        for i in orbital:
            line = [int(i) if i.isdigit() else i for i in i.strip().split()]
            if not line:
                break
            dats.append(line)
    return dats[1:]


def read_tdos(filename):
    dat = []
    with open(filename, "r") as tdos:
        for i in tdos:
            dat.append(
                np.asarray(i.strip().split(), dtype=float)
            )
    return np.asarray(dat)


def get_efermi(running_scf_log):
    for i in reverse_readfile(running_scf_log):
        if "EFERMI" in i:
            return round(
                float(i.strip().split('=')[-1].split()[0]), 4
            )


def get_bandgap(running_scf_log):
    for i in reverse_readfile(running_scf_log):
        if "Band Gap" in i:
            return round(
                float(i.strip().split('=')[-1].split()[0]), 4
            )


def get_energy(running_scf_log):
    for i in reverse_readfile(running_scf_log):
        if "!FINAL_ETOT_IS" in i:
            return round(
                float(
                    i.strip().split()[1]
                ), 4
            )


def get_total_magnetism(running_scf_log):
    for i in reverse_readfile(running_scf_log):
        if "total magnetism" in i:
            return round(
                float(
                    i.strip().split('=')[-1]
                ), 4
            )


def get_final_structure(out_dir):
    struct = sorted(
        [i for i in os.listdir(out_dir) if "STRU_ION" in i]
    )[-1]
    return os.path.join(out_dir, struct)


def read_band_dat(bands_x):
    bands = np.loadtxt(bands_x)
    _, k = bands.shape
    lb = np.asarray(list(range(1, k)))
    all_dat = []
    for item in bands:
        eig_val = np.vstack([lb, item[1:], [0, ] * (k - 1)])
        for j in eig_val.T:
            all_dat.append(j)
    return all_dat


class BandNscfLog:
    def __init__(self, filename):
        evf = False
        self._frac_kpts = []
        # self._eigenvalue = []
        self._kpoints = []
        self._nbands, self._nspin, self._nkstot, self._vbm, self._cbm, self._bandgap = [None, ] * 6
        with open(filename, "r") as nscf:
            for i in nscf:
                if "NBANDS" in i:
                    self._nbands = int(i.strip().split('=')[-1])
                if "SETUP K-POINTS" in i:
                    self._nspin = int(nscf.readline().strip('\n').split('=')[-1].strip())
                    self._nkstot = int(nscf.readline().strip('\n').split('=')[-1].strip())
                    for _ in range(2): nscf.readline()
                    for _ in range(self._nkstot):
                        self._frac_kpts.append(
                            np.asarray(nscf.readline().strip().split()[1:4], dtype=float)
                        )

                if "band eigenvalue in this processor" in i:
                    evf = True

                if evf:
                    if "k-points" in i:
                        kpt = np.asarray(i.strip().split(':')[-1].split(), dtype=float)
                        self._kpoints.append(kpt)
                    elif "final_state" in i:
                        # ev = np.asarray(i.strip().split()[1:], dtype=float)
                        # self._eigenvalue.append(ev)
                        continue
                    else:
                        if "Valence Band maximum is (eV):" in i:
                            self._vbm = float(i.strip().split('=')[-1])
                        if "Conduction Band minimum" in i:
                            self._cbm = float(i.strip().split('=')[-1])
                        if "Band Gap is" in i:
                            self._bandgap = float(i.strip().split('=')[-1])

    @property
    def nkpt(self):
        return self._nkstot

    @property
    def cbm(self):
        return self._cbm

    @property
    def vbm(self):
        return self._vbm

    @property
    def bandgap(self):
        return self._bandgap

    @property
    def nbands(self):
        return self._nbands

    @property
    def ispin(self):
        return self._nspin == 2

    def _get_cart_kpoints(self):
        if not self.ispin:
            return self._kpoints
        k = int(len(self._kpoints) / 2)

        return self._kpoints[:k]

    @property
    def frac_kpoints(self):
        return self._frac_kpts

    @property
    def cart_kpoints(self):
        return self._get_cart_kpoints()

    def get_eigenval(self, eigenvalue):
        n = self.nbands * self._nkstot
        if self.ispin:
            up, down = eigenvalue[:n], eigenvalue[n:]
            _eig = np.hstack([up, down])
        else:
            _eig = eigenvalue
        eig = [_eig[i: i + self.nbands] for i in range(0, n, self.nbands)]

        eigenval = _merge(eig, self.ispin, self._nkstot, self.nbands,
                          up_index=(1, 2), down_index=(4, 5))
        return eigenval


if __name__ == "__main__":
    pass
