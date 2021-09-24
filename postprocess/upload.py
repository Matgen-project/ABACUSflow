#!/usr/env/python3

from main import *
from pathlib import Path
import re
import os

from monty.os import cd
import pymongo

COMPATH=Path("/WORK/nscc-gz_material_1/ICSD_vasp/abacus_calc/matgen_scf/completed")


def yield_stru(root):
    for i in root.rglob('*'):
        if len(i.parts) == len(root.parts) + 1:
            yield comp(i), i

def comp(path):
    dos_out = path / "OUT.ABACUS"
    band_out = path / "BAND_OUT.ABACUS"
    scf_out = path / "SCF_OUT.ABACUS"
    stru = path / path.name
    if dos_out.exists() and band_out.exists() and scf_out.exists():
        files = {"stru": stru, "scf_log": scf_out / "running_scf.log",
                 "path": path, }
        return files
         

def get_stru(stru_filepath):
    return get_structure_dat(stru_filepath)

def get_mag(scf_dir):
    return get_magnetism(scf_dir)

def get_band(band_dir, stru_filename, scf_log_filepath):
    return get_bandstructure(band_dir, stru_filename, scf_log_filepath)


def get_dos(dos_dir):
    return get_density_of_states(dos_dir)

def get_paras(calc_dir):
    kpt = get_kpt(calc_dir, "scf")
    kpath = get_kpt(calc_dir, "band")
    return kpt, kpath

def get_dat(raw):
    db, stru_id, _ = re.split(r"[_|-]", raw["stru"].name)
    db += "_id"
    key = {db: int(stru_id)}
    stru = get_stru(raw["stru"])
    stru.update(key)
    band = get_band(raw["path"], raw["stru"].name, raw["scf_log"])
    band.update(key)
    dos = get_dos(raw["path"])
    dos.update(key)
    mag = get_mag(raw["path"])
    mag.update(key)
    cif = get_optimized_cif(raw["stru"])
    cif.update(key)
    return stru, mag, band, dos, cif, key

def get_db():
    addr = "12.11.70.140:10102"
    client = pymongo.MongoClient(addr)
    db = client["abacus_data"]
    return db

def upload_dat(db, *dat):
    stru, mag, band, dos, cif, key = dat
    stru_col = db["stru"] 
    band_col = db["bs_plotter"]
    dos_col = db["dos_plotter"]
    mag_col = db["mag"]
    cif_col = db["cif"]
    def _upload(col, data):
        exist = col.find_one(key)
        if exist is not None:
            col.update_one(exist, {'$set': data})
        else:
            col.insert_one(data)
    _upload(stru_col, stru)
    _upload(band_col, band)
    _upload(dos_col, dos)
    _upload(mag_col, mag)
    _upload(cif_col, cif)
    print(f"upload {key} sucessed.")

if __name__ == "__main__":
    adb = get_db()
    print(adb, " connected!")
    for calc, i in yield_stru(COMPATH):
        if calc is not None:
            try:
                res = get_dat(calc)
                upload_dat(adb, *res)
            except:
                os.system(f"echo {i} failed >> log") 
