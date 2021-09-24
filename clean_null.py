#!/bin/bash

import sys
from pathlib import Path
import os


def yield_null(file):
    with open(file, "r") as f:
        for line in f:
            if "NULL" in line:
                name = line.split()[0] 
                yield name



if __name__ == "__main__":
    args = sys.argv
    scf = Path("/WORK/nscc-gz_material_1/ICSD_vasp/abacus_calc/matgen_scf")
    r = args[1]
    for i in yield_null(r):
        stru = scf / i
        if stru.exists():
            os.system("rm -rf {}".format(stru))
