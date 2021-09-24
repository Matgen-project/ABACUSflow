#!/usr/bin/env python
# -*- coding: utf-8 -*-

from monty.os import cd

from core.structure import Crystal


class Inputs:
    def __init__(self, cur_dir):
        self.cur_dir = cur_dir

    def _read(self, filename, ftype):
        with cd(self.cur_dir):
            with open(filename, "r") as file:
                dat = file.read()

        return {ftype: dat}

    def get_KPT(self, ftype):
        _type = {"relax": "KPOINTS_relax", "scf": "KPOINTS_scf", "band": "KPATH"}
        if ftype == 'scf':
            name = "SCF_KPT"
        elif ftype == 'band':
            name = "BAND_KPT"
        else:
            name = "KPT"
        return self._read(name, _type.get(ftype))

    def get_structure(self, filename):
        with cd(self.cur_dir):
            return Crystal.matgen_structure_poscar_unopt(filename)

    def get_Input(self):
        return self._read("INPUT", "input")

    def get_potential(self, filename):
        with cd(self.cur_dir):
            cy = Crystal(filename).stru
            potential = dict(zip(cy.atom_species, cy.atom_potential))
            return {"potential_name": potential}


if __name__ == "__main__":
    pass
