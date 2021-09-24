#!/usr/env/python3

import sys
import os
import json
from pathlib import Path
from monty.os import cd
from subprocess import getoutput
from time import sleep

from ase.calculators.abacus.create_input import AbacusInput

HELPER = r"/WORK/nscc-gz_material_1/ICSD_vasp/abacus_calc/scripts/abacus_helper/"

def yield_stru(path):
   for i in path.rglob("*"):
       if len(i.parts) == len(path.parts) + 1:
           yield i

def _write_kpath(path):
    with cd(HELPER):
        kpath_dat = getoutput(f"python main.py -t kpath -s {path/path.name} -n 20")
    kpd = json.loads(kpath_dat)
    with cd(path):
        with open("KPT", "w") as f:
            for _, v in eval(kpath_dat).items():
                f.write(v)

def _write_input(path):
    with cd(path):
        nbands = getoutput("grep NBANDS SCF_OUT.ABACUS/running_scf.log |awk '{print $NF}' |head -n 1")
        ntype = getoutput("ls *.orb |wc -l")
        input = AbacusInput()
        input.set(atom_file=path.name,
                  ntype=int(ntype),
                  nbands=int(nbands),
                  kpoint_file='KPT',
                  pseudo_dir='./',
                  calculation='nscf',
                  nspin=2,
                  ecutwfc=60,
                  niter=50,
                  dr2=1.0e-9,
                  ethr=1.0e-7,
                  start_charge="file",
                  out_band=1,
                  smearing='gaussian',
                  sigma=0.02,
                  ks_solver='genelpa',
                  basis_type='lcao',
                  mixing_type='pulay',
                  mixing_beta=0.4,
                  gamma_only=0)
        input.write_input_input()
        os.system('yhbatch -N 1 abacus.sh') 
        sleep(1)


def make_band_inputs(scf):
    with cd(scf):
        os.system("mv *.out ./OUT.ABACUS/")
        os.system("mv OUT.ABACUS SCF_OUT.ABACUS")
        os.system("mv INPUT SCF_INPUT")
        os.system("mv KPT SCF_KPT")
    _write_kpath(scf)
    _write_input(scf)


if __name__ == "__main__":
    args = sys.argv

    scf = Path(args[1])
    for s in yield_stru(scf):
        s = s.absolute()
        try:
            make_band_inputs(s)
        except:
            os.system(f"echo {s} >> err.log")    
   
