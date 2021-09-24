#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os
from functools import reduce

from utils.output import read_tdos, read_pdos, read_orbital

# np.set_printoptions(precision=4)


class DensityOfStates:
    def __init__(self, dos_dir):
        out_dir = os.path.join(dos_dir, "OUT.ABACUS")
        self.out = out_dir
        self._fp = os.path.join(self.out, "PDOS")
        self._orb = os.path.join(self.out, "Orbital")
        self.ispin = self._get_spin()

    def _get_spin(self):
        with open(self._fp, "r") as f:
            for i in f:
                if "nspin" in i:
                    spin = int(i.strip().replace("<nspin>", '').replace("</nspin>", ''))
        if spin == 1:
            return False
        return True

    def _get_tdos(self):
        dos1_smearing = os.path.join(self.out, "DOS1_smearing.dat")
        if self.ispin:
            self.tdos_label = ["energy", "dos_up", "dos_down", "int_dos_up", "int_dos_down"]
            dos2_smearing = os.path.join(self.out, "DOS2_smearing.dat")
            up, down = read_tdos(dos1_smearing), read_tdos(dos2_smearing)
            energy = up[:, 0]
            up_dos_dat = up[:, 1]
            up_int_dos_dat = up[:, 2]
            down_dos_dat = down[:, 1] * -1
            down_int_dos_dat = down[:, 2]
            tdos = np.vstack([energy, up_dos_dat, down_dos_dat, up_int_dos_dat, down_int_dos_dat])
            return tdos.T
        else:
            self.tdos_label = ["energy", "dos", "int_dos"]
            tdos_data = read_tdos(dos1_smearing)
            return tdos_data

    def _get_matgen_tdos_old_style(self):
        tdos_dat = self._get_tdos()
        return {
            "TDOS": {
                "tdos_labels": self.tdos_label,
                "tdos_data": tdos_dat.tolist()}

        }

    def _get_orb(self):
        _orb = read_orbital(self._orb)
        info = {}
        for i in _orb:
            io, spec, l, m, z, sym = i
            if info.get(spec) is None:
                info[spec] = [(l, m, z, sym)]
            else:
                _info = info.get(spec)
                _info.append((l, m, z, sym))
                info[spec] = _info

        return info

    def _get_pdos(self):
        dat, energy = read_pdos(self._fp)
        orbs = self._get_orb()
        pdos_dict = {}
        for item in dat:
            spec = item["species"]
            ai = item["atom_index"]
            name = spec + "_" + str(ai)
            l, m, z = list(map(item.get, ['l', 'm', 'z']))
            spec_orb = orbs[spec]
            pdos = item["pdos"]
            for i in spec_orb:
                if (l, m, z) == i[:3]:
                    sym = i[3]
                    if pdos_dict.get(name) is None:
                        pdos_dict[name] = [[(l, m, z, sym), pdos]]
                    else:
                        pdos_dict[name].append([(l, m, z, sym), pdos])

        return pdos_dict, energy

    def _get_matgen_pdos_old_style(self):
        tmp, energy = self._get_pdos()
        pd = {}
        for k, v in tmp.items():
            m = {}
            for item in v:
                sym, pdos = item
                lb = sym[-1][0]
                if m.get(lb) is None:
                    m[lb] = pdos
                else:
                    m[lb] += pdos
            pd[k] = m
        _plus_pd = {}
        for m, n in pd.items():
            spec = m.split('_')[0]
            if _plus_pd.get(spec) is None:
                _plus_pd[spec] = n
            else:
                for orb, orb_pd in n.items():
                    _plus_pd[spec][orb] += orb_pd

        def add(x, y):
            return x + y

        def _vstack(x, y):
            return np.vstack([x, y])

        pdos_dict = {}
        for element, pd in _plus_pd.items():
            up_lst, down_list = [], []
            pdos_label = list(pd.keys())
            pdos_label.append("total")
            pdos_label.insert(0, "Energy(eV)")
            if self.ispin:
                up_name = element + "_up"
                down_name = element + "_down"
                for orb, ipd in pd.items():
                    up, down = ipd[:, 0], ipd[:, 1]
                    up_lst.append(up)
                    down_list.append(down * -1)
                total_up = reduce(add, up_lst)
                total_down = reduce(add, down_list)
                up_lst.append(total_up)
                up_lst.insert(0, energy)
                down_list.append(total_down)
                down_list.insert(0, energy)
                up_data = reduce(_vstack, up_lst).T
                down_data = reduce(_vstack, down_list).T
                up_dict = {up_name: {"pdos_label": pdos_label, "pdos_data": up_data.tolist()}}
                down_dict = {down_name: {"pdos_label": pdos_label, "pdos_data": down_data.tolist()}}
                pdos_dict.update(up_dict)
                pdos_dict.update(down_dict)
            else:
                for orb, ipd in pd.items():
                    up_lst.append(ipd[:, 0])
                total_up = reduce(add, up_lst)
                up_lst.append(total_up)
                up_lst.insert(0, energy)
                up_data = reduce(_vstack, up_lst).T
                up_dict = {element: {"pdos_label": pdos_label, "pdos_data": up_data.tolist()}}
                pdos_dict.update(up_dict)

        return {"PDOS": pdos_dict}

    def get_dos_using_matgen_old_fmt(self):
        tdos = self._get_matgen_tdos_old_style()
        pdos = self._get_matgen_pdos_old_style()

        tdos.update(pdos)
        return {"Density_of_states": tdos}

    def plot(self):
        import matplotlib.pyplot as plt
        dos_data = self.get_dos_using_matgen_old_fmt()
        pdos_data = dos_data.get('Density_of_states').get('PDOS')
        all_atoms = list(pdos_data.keys())
        for at in all_atoms:
            plt.xlabel("Energy(eV)")
            plt.ylabel("electrons/eV")
            plt.title(at.split('_')[0] + "_Density_of_states")
            dlb = pdos_data.get(at).get('pdos_label')
            dd = pdos_data.get(at).get('pdos_data')
            Energy = dd.T[0]
            for i in range(1, len(dlb)):
                plt.plot(Energy, dd.T[i])
            if at.split('_')[-1] == 'up':
                continue
            plt.show()
        plt.title("Total_Density_of_states")
        plt.xlabel("Energy(eV)")
        plt.ylabel("electrons/eV")
        tdos_data = dos_data.get('Density_of_states').get('TDOS')
        tlb = tdos_data.get('tdos_labels')
        td = tdos_data.get('tdos_data')
        if len(tlb) > 3:
            Energy, dos_up, dos_down = td.T[:3]
            plt.plot(Energy, dos_up, color='red')
            plt.plot(Energy, dos_down, color='blue')
        else:
            Energy, dos_up = td.T[:2]
            plt.plot(Energy, dos_up, color='red')
        plt.show()


if __name__ == "__main__":
    pass
