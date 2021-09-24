#!/usr/env/python3

import json

with open(r"abacus_simple.json", "r") as f:
    si = json.load(f)


def support():
    pfile = r"../../abacus.pot"
    ofile = r"../../abacus.orb"
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

alls = support()
for k in alls:
    if k not in list(si.keys()):
        print(k)
     
