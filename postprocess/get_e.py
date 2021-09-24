#!/usr/env/python3

from main import *
from pathlib import Path
import re
import os
from subprocess import getoutput

from monty.os import cd
import pymongo

COMPATH=Path("/WORK/nscc-gz_material_1/ICSD_vasp/abacus_calc/matgen_scf/completed")
ERRORS = Path("/WORK/nscc-gz_material_1/ICSD_vasp/abacus_calc/matgen_scf/some_errors")


def yield_stru(root):
    for i in root.rglob('*'):
        if len(i.parts) == len(root.parts) + 1:
            yield comp(i), i

def comp(path):
    scf_out = path / "SCF_OUT.ABACUS"
    stru = path / path.name
    if scf_out.exists():
        files = {"stru": stru, "scf_log": scf_out / "running_scf.log",
                 "path": path, "scf": scf_out }
        return files

def get_energy(log, stru, mid):
    print(log)
    is_sim = False
    cmd1 = f"grep \'!FINAL_ETOT_IS\' {log}"
    total_e = float(getoutput(cmd1).split()[1])
    cmd2 = f"grep \'TOTAL ATOM NUMBER\' {log}"
    n = int(getoutput(cmd2).split()[-1])
    cmd3 = f"grep \'ntype\' {log}"
    ntype = int(getoutput(cmd3).split()[-1])
    cmd4 = f"grep EFERMI {log}"
    ef = getoutput(cmd4).split()[2]
    if ntype == 1:
        is_sim = True
    total_e_pa = total_e / n
    cmd5 = f"grep \'atom label for species\' {log}"
    syb = []
    lbs = getoutput(cmd5)
    for line in lbs.split("\n"):
        syb.append(line.split("=")[-1].strip(' '))
    cmd6 = f"grep \'number of atom for this type\' {log}"
    ns = getoutput(cmd6)
    num = []
    for yl in ns.split('\n'):
        num.append(int(yl.split("=")[-1]))
    species = dict(zip(syb, num))

    return {"id":mid, "energy": total_e, "epa": total_e_pa, "efermi": ef,"is_sim": is_sim, "symbol": species} 
         

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
    #stru = get_stru(raw["stru"])
    #stru.update(key)
    #band = get_band(raw["path"], raw["stru"].name, raw["scf_log"])
    #band.update(key)
    #dos = get_dos(raw["path"])
    #dos.update(key)
    #mag = get_mag(raw["path"])
    #mag.update(key)
    #cif = get_optimized_cif(raw["stru"])
    #cif.update(key)
    #return stru, mag, band, dos, cif, key
    e = get_energy(raw["scf_log"], raw["stru"], int(stru_id))
    e.update(key)
    return e, key

def goe(e, fe, k):
    r = {}
    r.update(k)
    r['efermi'] = e['efermi']
    r['energy'] = e['energy']
    r.update(fe)
    return r


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

def upload_eng(db, dat):
    en_col = db["energy"] 
    def _upload(col, data):
        exist = col.find_one(key)
        if exist is not None:
            col.update_one(exist, {'$set': data})
        else:
            col.insert_one(data)
    _upload(en_col, dat)
    print(f"upload sucessed.")


def calcfe(item, s):
    if item['is_sim'] == 'True':
        syb = list(item['symbol'].keys())[0]
        te = float(item['energy'])
        s2_id = s[syb][1]

        if int(item["id"]) == int(s2_id):
            foe = 0
        else:
            s2_e = s[syb][0]
            v = int(list(item['symbol'].values())[0])
            foe = (te - v * s2_e) / v
    else:
        syb = item['symbol']
        te = float(item['energy'])
        an = 0
        for a, v in syb.items():
            ie = int(v) * s[a][0]
            te -= ie
            an += int(v)
        foe = round(te / an, 4)

    #return f"{item['id']}\t{foe}"
    return {"formation_energy": foe}


if __name__ == "__main__":
    import json
    adb = get_db()
    print(adb, " connected!")
    c = 0
    with open("abacus_simple.json", "r") as f:
        s = json.load(f)
    for calc, i in yield_stru(COMPATH):
        if calc is not None:
            try:
                res, key = get_dat(calc)
                foe = calcfe(res, s)
            except Exception as e:
                c += 1
                #os.system(f"mv {i} {ERRORS}")
                continue
            else:
                ans = goe(res, foe, key)
                print(ans)
                upload_eng(adb, ans)
    print("errors: ", c)
            
     
                        
