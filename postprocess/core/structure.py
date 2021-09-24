#!/usr/bin/env python
# -*- coding: utf-8 -*-


from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io import cif

from utils.stru import Stru
from utils.kpt import Kpt


class Crystal:
    def __init__(self, filename_stru):
        self.fn = filename_stru
        self.stru = Stru.from_stru(filename_stru)
        self.crystal = self.convert()

    def _to_poscar(self):
        e_str = "  ".join(self.stru.uniq_atom_species)
        an_str = "  ".join([str(i) for i in self.stru.atom_number])
        pos_str = ""
        for i in self.stru.frac_coords:
            pos_str += "   ".join(str(j) for j in i)
            pos_str += "\n"
        formula = self.stru.formula
        la, lb, lc = self.stru.lattice.lattice
        las = "   ".join(str(i) for i in la.tolist())
        lbs = "   ".join(str(j) for j in lb.tolist())
        lcs = "   ".join(str(k) for k in lc.tolist())
        poscar_fmt = f"{formula}\n" \
                     f"1.0\n" \
                     f" {las}\n" \
                     f" {lbs}\n" \
                     f" {lcs}\n" \
                     f" {e_str}\n" \
                     f" {an_str}\n" \
                     f"Direct\n" \
                     f"{pos_str}"

        return poscar_fmt

    def convert(self):
        return Structure.from_str(self._to_poscar(), fmt="poscar")

    def get_kpath(self, n=20):
        kpath = Kpt.generate_kpath_from(self.fn)
        point_coords = kpath['point_coords']
        _path = kpath['path']
        plst = []
        count = 0
        for i, ipath in enumerate(_path):
            if i > 0:
                last_k = _path[i - 1][1]
                if last_k == ipath[0]:
                    ipath = [ipath[1]]
            for _p in ipath:
                coord = point_coords.get(_p)
                line = "%.4f  %.4f  %.4f  %d  #  %s\n" % (coord[0], coord[1], coord[2], n, _p)
                plst.append(line)
                count += 1
        _f = plst[-1].strip().split()
        _f[3] = '1'
        _fs = "  ".join(_f)
        pstring = ''.join(plst[:-1]) + _fs
        head = f"K_POINTS \n" \
               f"{count} # number of high symmetry lines\n" \
               "Line # line-mode\n"
        return {"kpath": head + pstring}

    def matgen_structure_old_style(self):
        csga = SpacegroupAnalyzer(self.crystal)
        conv_cell = csga.get_conventional_standard_structure()
        formula = conv_cell.composition.formula
        ccl = conv_cell.as_dict().get('lattice')
        ccspg = conv_cell.get_space_group_info()
        prmi_cell = csga.get_primitive_standard_structure()
        pcl = prmi_cell.as_dict().get('lattice')
        cs = csga.get_crystal_system()
        pgs = csga.get_point_group_symbol()
        c_cell_elements = conv_cell.sites
        p_cell = conv_cell.get_primitive_structure()
        p_cell_elements = p_cell.sites
        density = self.crystal.density

        def get_sp_and_coor(elements):
            sp, coor = [], []
            for i in elements:
                coor.append(i.frac_coords.tolist())
                sp.append(i.species_string)
            return {'atoms_order': sp, 'atoms_coordinates': coor}

        c_cell_res = get_sp_and_coor(c_cell_elements)
        p_cell_res = get_sp_and_coor(p_cell_elements)

        return {"formula": formula, "conventional_cell": ccl, "conventional_cell_site": c_cell_res,
                "primitive_cell": pcl, "primitive_cell_site": p_cell_res, "crystal_system": cs,
                "density": density, "point_group": pgs, "spacegroup": ccspg,
                "elements": self.stru.atom_species}

    def matgen_structure_cif_opt(self):

        cf_data = cif.CifWriter(self.crystal)
        cf = str(cf_data).replace('# generated using pymatgen',
                                  '# geometry optimization by matgen')

        return {'cif_data': cf}

    @classmethod
    def matgen_structure_poscar_unopt(cls, filename):
        unopt_structure = cls(filename)._to_poscar()
        return {"POSCAR": unopt_structure}


if __name__ == "__main__":
    pass
