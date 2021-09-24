#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os

from utils.stru import Stru
from utils.kpt import Kpt
from utils.output import BandNscfLog, get_efermi, read_band_dat
from utils._base.spin import Spin

# np.set_printoptions(precision=5)


class BandStructureSymmLine:
    def __init__(self, band_dir, stru_filename, scf_log):
        stru_filepath = os.path.join(band_dir, stru_filename)
        #self.out_dir = os.path.join(band_dir, "OUT.ABACUS")
        self.out_dir = os.path.join(band_dir, "BAND_OUT.ABACUS")
        self.lattice = Stru.from_stru(stru_filepath).lattice.reciprocal_lattice()
        self._b1 = os.path.join(self.out_dir, "BANDS_1.dat")
        #kpt_filepath = os.path.join(band_dir, "KPT")
        kpt_filepath = os.path.join(band_dir, "BAND_KPT")
        coords, self.ikpt, lbs, *_ = Kpt.from_kpt(kpt_filepath)
        self.kpoints = []
        labels_dict = dict(zip(lbs, coords))
        nscf_log = os.path.join(self.out_dir, "running_nscf.log")
        self.log = BandNscfLog(nscf_log)
        self.scf_log = scf_log
        _kpt = self.log.frac_kpoints
        for k in _kpt:
            # let see if this kpoint has been assigned a label
            label = None
            for c in labels_dict:
                c_val = labels_dict[c]
                if isinstance(c_val, Kpt):
                    c_val = c_val.frac_coords
                if np.linalg.norm(k - c_val) < 0.0001:
                    label = c
                    labels_dict[label] = Kpt(
                        k,
                        self.lattice,
                        label=label,
                        coords_are_cartesian=False,
                    )
            self.kpoints.append(
                Kpt(
                    k, self.lattice, label=label, coords_are_cartesian=False
                )
            )

        self.high_kpts, self.distance = self._get_hk()

    def _get_hk(self):
        distance = []
        tick_distance = []
        tick_labels = []
        previous_kpoint = self.kpoints[0]
        previous_distance = 0.0
        previous_label = self.kpoints[0].label
        for i in range(len(self.kpoints)):
            label = self.kpoints[i].label
            if label is not None and previous_label is not None:
                distance.append(previous_distance)
            else:
                distance.append(
                    np.linalg.norm(
                        self.kpoints[i].cart_coords - previous_kpoint.cart_coords
                    )
                    + previous_distance
                )
            previous_kpoint = self.kpoints[i]
            previous_distance = distance[i]
            previous_label = label
            if label:
                tick_distance.append(distance[i])
                tick_labels.append(label)

        high_kpts = {'distance': tick_distance, 'label': tick_labels}

        return high_kpts, distance

    @property
    def bandgap(self):
        return self.log.bandgap

    @property
    def cbm(self):
        return self.log.cbm

    @property
    def vbm(self):
        return self.log.vbm

    @property
    def ispin(self):
        return self.log.ispin

    @property
    def nbands(self):
        return self.log.nbands

    def _read_eigenval(self):
        b1d = read_band_dat(self._b1)
        if self.ispin:
            self._b2 = os.path.join(self.out_dir, "BANDS_2.dat")
            b2d = read_band_dat(self._b2)
            b1d.extend(b2d)

        return np.asarray(b1d)

    @property
    def bands(self):
        return self.log.get_eigenval(self._read_eigenval())

    @property
    def nkpt(self):
        return self.log.nkpt

    @property
    def efermi(self):
        return get_efermi(self.scf_log)

    def is_metal(self, efermi_tol=1e-4):
        """
        Check if the band structure indicates a metal by looking if the fermi
        level crosses a band.

        Returns:
            True if a metal, False if not
        """
        band = []
        for i in range(self.nbands):
            for item in self.bands[Spin.up]:
                band.append(item[i][0])
            tmp = np.asarray(band)
            band = []
            if np.any(tmp < -efermi_tol) and np.any(
                    tmp > efermi_tol
            ):
                return True
        return False

    def get_bandstructure_using_matgen_old_fmt(self):
        energy_data = {}
        if self.ispin:
            labels = ["Wave_vector", "spin_up", "spin_down"]
            up, down = self.bands[Spin.up], self.bands[Spin.down]
            for band_index in range(self.nbands):
                key = f"band_index_{band_index + 1}"
                vals = []
                for i in range(self.nkpt):
                    up_band_dat = up[i][band_index]
                    down_band_dat = down[i][band_index]
                    vals.append([self.distance[i], up_band_dat[0], down_band_dat[0]])
                energy_data.update({key: vals})
        else:
            labels = ["Wave_vector", "Energy_level"]
            up = self.bands[Spin.up]
            for band_index in range(self.nbands):
                key = f"band_index_{band_index + 1}"
                vals = []
                for i in range(self.nkpt):
                    up_band_dat = up[i][band_index]
                    vals.append([self.distance[i], up_band_dat[0]])
                energy_data.update({key: vals})
        kpath = {'High_Kpoints_labels': self.high_kpts['label'],
                 'High_Kpoints_coordinates': self.high_kpts['distance']}

        return {'Band_Structure': {
            'Band_Gap': self.bandgap if not self.is_metal() else 0.0,
            'Energy_data_labels': labels,
            'Energy_data': energy_data,
            'Spin_state': self.ispin,
            'Hk_points': kpath}}

    def plot(self):
        import matplotlib.pyplot as plt
        band_data = self.get_bandstructure_using_matgen_old_fmt()
        bd_lbs = band_data.get('Band_Structure').get('Energy_data_labels')
        bd_arrary = band_data.get('Band_Structure').get('Energy_data')
        hklbs = band_data.get('Band_Structure').get('Hk_points')
        lbcoord = hklbs.get('High_Kpoints_coordinates')
        plt.title("Band_Structure")
        plt.xlabel("Wave_vector")
        plt.ylabel("Energy(eV)")
        plt.ylim(-15, 15)
        for bandi, bandd in bd_arrary.items():
            if len(bd_lbs) == 3:
                x, y1, y2 = np.asarray(bandd).T[0], np.asarray(bandd).T[1], np.asarray(bandd).T[2]
                plt.plot(x, y1, color='black')
                plt.plot(x, y2, color='blue')
            else:
                x, y = np.asarray(bandd).T[0], np.asarray(bandd).T[1]
                plt.plot(x, y, color='blue')
        for k, i in enumerate(lbcoord):
            plt.axvline(i, color='red')
        plt.show()
        # print(self.is_metal())


if __name__ == "__main__":
    pass

