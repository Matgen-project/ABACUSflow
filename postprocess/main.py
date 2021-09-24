#!/usr/bin/env python
# -*- coding: utf-8 -*-


from core.bandstructure import BandStructureSymmLine
from core.dos import DensityOfStates
from core.structure import Crystal
from core.optimize import CellRelax, Scf
from core.inputs import Inputs


def get_bandstructure(band_dir, stru_filename, scf_log_filepath):
    band = BandStructureSymmLine(band_dir, stru_filename, scf_log_filepath)
    band_data = band.get_bandstructure_using_matgen_old_fmt()
    # print(band_data)
    # band.plot()
    return band_data


def get_density_of_states(dos_dir):
    dos = DensityOfStates(dos_dir)
    dos_data = dos.get_dos_using_matgen_old_fmt()
    # print(dos_data)
    # dos.plot()
    return dos_data


def get_magnetism(scf_dir):
    mag_data = Scf(scf_dir).get_mag()
    return mag_data


def get_optimized_stru_filepath(relax_dir):
    return CellRelax(relax_dir).get_result()


def get_efermi(scf_dir):
    return Scf(scf_dir).get_efermi()


def get_energy(scf_dir):
    return Scf(scf_dir).get_energy()


def get_optimized_cif(stru_filepath):
    return Crystal(stru_filepath).matgen_structure_cif_opt()


def get_unoptimized_poscar(stru_filepath):
    return Crystal.matgen_structure_poscar_unopt(stru_filepath)


def get_structure_dat(stru_filepath):
    return Crystal(stru_filepath).matgen_structure_old_style()


def get_kpath(stru_filepath, n=20):
    return Crystal(stru_filepath).get_kpath(n)


def get_kpt(calc_dir, ktype):
    return Inputs(calc_dir).get_KPT(ktype)


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='Abacus Post-process tool',
    )
    parser.add_argument('-t', '--type')
    parser.add_argument('-d', '--dir')
    parser.add_argument('-s', '--stru')
    parser.add_argument('-l', '--log')
    parser.add_argument('-n', '--num', type=int, default=20)
    parser.add_argument('-k', '--ktype', choices=['relax', 'scf', 'band'])
    args = parser.parse_args()
    if args.type == 'sp':
        res = get_structure_dat(args.stru)
    elif args.type == 'band':
        res = get_bandstructure(args.dir, args.stru, args.log)
    elif args.type == 'dos':
        res = get_density_of_states(args.dir)
    elif args.type == 'mag':
        res = get_magnetism(args.dir)
    elif args.type == 'energy':
        res = get_energy(args.dir)
    elif args.type == 'efermi':
        res = get_efermi(args.dir)
    elif args.type == 'cif':
        res = get_optimized_cif(args.stru)
    elif args.type == 'poscar':
        res = get_unoptimized_poscar(args.stru)
    elif args.type == 'opt':
        res = get_optimized_stru_filepath(args.dir)
    elif args.type == 'kpath':
        res = get_kpath(args.stru, args.num)
    elif args.type == 'kpts':
        res = get_kpt(args.dir, args.ktype)
    else:
        raise TypeError

    #print(res)
    return res


if __name__ == "__main__":
    """
        # BandStructure
    stru = "icsd_9852-Ti2O4"
    band_spin1 = r"./test/icsd_9852-Ti2O4/SPIN1"
    spin1_scf_log = r"./test/icsd_9852-Ti2O4/SPIN1/OUT.ABACUS/running_scf.log"
    band_spin2 = r"./test/icsd_9852-Ti2O4/SPIN2"
    spin2_scf_log = r"./test/icsd_9852-Ti2O4/SPIN2/OUT.ABACUS/running_scf.log"
    # x = get_bandstructure(band_spin1, stru, spin1_scf_log)
    # y = get_bandstructure(band_spin2, stru, spin2_scf_log)
    # DOS
    # spin1
    spin1 = r"./test/icsd_23076-Sr1Ti1O3/DOS_SPIN1/"
    # s1 = get_density_of_states(spin1)
    spin2 = r"./test/icsd_23076-Sr1Ti1O3/DOS_SPIN2/"
    # s2 = get_density_of_states(spin2)
    # Magnetism
    sd = r"./test/icsd_23076-Sr1Ti1O3/SCF"
    m = get_magnetism(sd)
    print(m)
    e = get_energy(sd)
    print(e)
    stru_fp = band_spin1 + '/' + stru
    cif = get_optimized_cif(stru_fp)
    print(cif)
    un_opt_poscar = get_unoptimized_poscar(stru_fp)
    print(un_opt_poscar)
    sp = get_structure_dat(stru_fp)
    print(sp)
    """
    main()
