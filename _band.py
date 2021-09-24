#!/bin/env/python
wl = []
with open(r"./results.dat", "r") as f:
    for dat in f.readlines():
        dat_line = dat.split()
        want = dat_line[0] + ' ' + dat_line[-1]
        wl.append(want)

print(wl)

with open(r"./abacus.bandgap", "w") as abg:
    for item in wl[1:]:
        abg.write(item + '\n')
