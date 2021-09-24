#!/usr/env/python3


from pathlib import Path
import shutil

def support():
    pfile = r"./abacus.pot"
    ofile = r"./abacus.orb"
    def read(file):
        su = set()
        with open(file, "r") as s:
            for i in s:
                su.add(i.strip("\n"))
        return su
    psu = read(pfile)
    osu = read(ofile)
    print("Pot support: ", len(psu))
    print("Orb support: ", len(osu))
    return psu & osu

def yield_stru():
    file = Path(r"/WORK/nscc-gz_material_1/ICSD_vasp/abacus_calc/matgen_scf/simple_substance/stru_from_matgen")
    for stru in file.rglob("*.vasp"):
        yield stru


def get_e(stru):
    with open(stru, "r") as f:
        for idx, line in enumerate(f):
            if idx == 5:
                elements = line.strip('\n').split()
    return elements

def in_support(es, sue):
    for i in es:
        if i not in sue:
            return False
    return True


if __name__ == "__main__":
    sue = support()
    print(sue)
    print("total support: ", len(sue))
    sim = set()
    for stru in yield_stru():
        ye = get_e(stru)
        print(ye) 
        if len(ye) == 1:
            sim.add(ye[0])
        if in_support(ye, sue):
            shutil.copy(stru, r"/WORK/nscc-gz_material_1/ICSD_vasp/abacus_calc/matgen_scf/simple_substance/run")
    print(sim)
    print(len(sim))
            




