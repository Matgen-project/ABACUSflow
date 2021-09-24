#!~/software/anaconda3/bin/python

import os
import sys

args =sys.argv

STRUPATH=args[1]
CALC_PATH = args[2]

with open('results.dat', 'w') as f:
    f.write('{0:<25s} {1:<25s} {2:<25s} {3:<25s}\n'.format(
            '# StructureName',
            'Energy(eV)',
            'Magnetism(Bohr mag/cell)',
            'Band gap(eV)'))


for root, dirs, files in os.walk(top=STRUPATH, topdown=True):
    for file in files:
        nm = file.split('.')[0]
        name = os.path.join(CALC_PATH, nm)
        print(name)
        #cmd1 = ("grep '!FINAL_ETOT_IS' " + name + "/OUT.ABACUS/running_scf.log | awk '{print $2} ' > tmpE") 
        cmd1 = ("grep '!FINAL_ETOT_IS' " + name + "/SCF_OUT.ABACUS/running_scf.log | awk '{print $2} ' > tmpE") 
        cmd2 = ("grep 'total magnetism' " + name + 
                "/OUT.ABACUS/running_scf.log | tail -n 1 | awk '{print $6}' > tmpMag")
        cmd3 = ("grep 'Band Gap' " + name + "/OUT.ABACUS/running_scf.log | awk '{print $6} ' > tmpBandGap") 
        os.system(cmd1)
        os.system(cmd2)
        os.system(cmd3)
        with open('tmpE', 'r') as f:
            ele = f.readline()
            if(ele == ''):
                energy = 'NULL'
            else:
                energy = ele.split()[0]
        with open('tmpMag', 'r') as f:
            ele = f.readline()
            if (ele == ''):
                mag = 'NULL'
            else:
                mag = ele.split()[0]
        with open('tmpBandGap', 'r') as f:
            ele = f.readline()
            if (ele == ''):
                gap = 'NULL'
            else:
                gap = ele.split()[0]

        with open('results.dat', 'a') as f:
            f.write('{0:<25s} {1:<25s} {2:<25s} {3:<25s}\n'.format(nm, energy, mag, gap))
