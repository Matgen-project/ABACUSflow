#!/usr/bin/env python
# -*- coding: utf-8 -*-

from os.path import *
from utils.output import (
    get_energy, get_efermi, get_total_magnetism,
    get_bandgap, get_final_structure
)


class CellRelax:
    def __init__(self, cell_relax_dir):
        self._relax_dir = join(cell_relax_dir, "OUT.ABACUS")

    def get_result(self) -> str:
        return get_final_structure(self._relax_dir)


class Relax(CellRelax):
    def __init__(self, relax_dir):
        super(Relax, self).__init__(relax_dir)


class Scf:
    def __init__(self, scf_dir):
        self._scf_log = join(scf_dir, "OUT.ABACUS/running_scf.log")

    def get_efermi(self):
        return {"efermi": get_efermi(self._scf_log)}

    def get_energy(self):
        return {"final": get_energy(self._scf_log)}

    def get_mag(self):
        return {"total": get_total_magnetism(self._scf_log)}

    def get_bandgap(self):
        return get_bandgap(self._scf_log)

    @property
    def log(self):
        return self._scf_log


if __name__ == "__main__":
    pass

