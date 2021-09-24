# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 10:31:30 2018

@author: shenzx
"""
from os.path import join, basename, exists
from ase import Atoms
import numpy as np
import shutil
from ase.calculators.abacus.potential import PotentialDict
from ase.calculators.abacus.basis import BasisDict
# from abacus.potential import PotentialDict
# from abacus.basis import BasisDict
# import sys
# sys.path.append("E:\Git\project\ase-abacus\ase-abacus")


def potential_list():
    return list(PotentialDict.keys())


def basis_list():
    return list(BasisDict.keys())


def judge_exist_stru(stru=None):
    if stru is None:
        return False
    else:
        return True


def read_ase_stru(stru=None, coordinates_type="Cartesian"):
    if judge_exist_stru(stru):
        atoms_list = []
        atoms_position = []
        atoms_masses = []
        atoms_magnetism = []
        atoms_all = stru.get_chemical_symbols()

        # sort atoms according to atoms
        for atoms_all_name in atoms_all:
            temp = True
            for atoms_list_name in atoms_list:
                if atoms_all_name == atoms_list_name:
                    temp = False
                    break

            if temp:
                atoms_list.append(atoms_all_name)

        for atoms_list_name in atoms_list:
            atoms_position.append([])
            atoms_masses.append([])
            atoms_magnetism.append([])

        # get position, masses, magnetism from ase atoms
        if coordinates_type == 'Cartesian':
            for i in range(len(atoms_list)):
                for j in range(len(atoms_all)):
                    if atoms_all[j] == atoms_list[i]:
                        atoms_position[i].append(list(
                                stru.get_positions()[j]))
                        atoms_masses[i] = stru.get_masses()[j]
                        atoms_magnetism[i] = stru.get_initial_magnetic_moments()[j]
                        # update 20201230

        elif coordinates_type == 'Direct':
            for i in range(len(atoms_list)):
                for j in range(len(atoms_all)):
                    if atoms_all[j] == atoms_list[i]:
                        atoms_position[i].append(list(
                                stru.get_scaled_positions()[j]))
                        atoms_masses[i] = stru.get_masses()[j]
                        atoms_magnetism[i] = stru.get_initial_magnetic_moments()[j]
                        # update 20201230
        else:
            raise ValueError("'coordinates_type' is ERROR,"
                             "please set to 'Cartesian' or 'Direct'")

        return atoms_list, atoms_masses, atoms_position, atoms_magnetism


def set_potential(atoms_list=None, pseudo_dir="./", potential_name=None):
    if atoms_list is None:
        print(" Please set right 'atoms_list' ")

    else:
        potential = []
        if potential_name is None:
            for atoms_list_name in atoms_list:
                potential.append(
                        join(pseudo_dir,
                             PotentialDict['PotLDA'][atoms_list_name]))

        elif type(potential_name) == str:
            PotList = potential_list()
            if potential_name in PotList:
                for atoms_list_name in atoms_list:
                    potential.append(join(
                            pseudo_dir,
                            PotentialDict[potential_name][atoms_list_name]))
            else:
                raise ValueError("'potential_name' is ERROR")
            """
            if potential_name == 'PotLDA':
                for atoms_list_name in atoms_list:
                    potential.append(
                            join(pseudo_dir,
                                 PotentialDict['PotLDA'][atoms_list_name]))

            elif potential_name == 'PotPBE':
                for atoms_list_name in atoms_list:
                    potential.append(join(
                            pseudo_dir,
                            potential.PotPBE[atoms_list_name]))

            elif potential_name == 'PotPBESG15':
                for atoms_list_name in atoms_list:
                    potential.append(join(
                            pseudo_dir,
                            potential.PotPBESG15[atoms_list_name]))

            else:
                raise ValueError("'potential_name' is ERROR")
        """

        elif type(potential_name) == list:
            ele_name = {}
            for i in potential_name:
                with open(join(pseudo_dir, i), 'r') as f:
                    lines = f.readlines()

                for line in lines:
                    line = line.replace('=', ' = ')
                    line = line.replace('"', ' " ')
                    data = line.split()
                    if len(data) == 0:
                        continue

                    elif data[0] == 'element':
                        ele_name[data[3]] = i
                        break

                    elif len(data) == 2 and data[1] == 'Element':
                        ele_name[data[0]] = i
                        break

                    else:
                        continue

            for atoms_list_name in atoms_list:
                potential.append(join(pseudo_dir,
                                      ele_name[atoms_list_name]))

        else:
            raise ValueError("Please sure what you do!!! ")

        return potential


def set_basis(atoms_list=None, basis_dir="./", basis_name=None):
    if atoms_list is None:
        print(" Please set right 'atoms_list'")

    else:
        basis = []
        if basis_name is None:
            for atoms_list_name in atoms_list:
                basis.append(join(
                        basis_dir,
                        BasisDict['LDAmin'][atoms_list_name]))

        elif type(basis_name) == str:
            BasisList = basis_list()
            if basis_name in BasisList:
                for atoms_list_name in atoms_list:
                    basis.append(join(
                            basis_dir,
                            BasisDict[basis_name][atoms_list_name]))
            else:
                raise ValueError("'basis_name' is ERROR")

        elif type(basis_name) == list:
            ele_name = {}
            for i in basis_name:
                with open(join(basis_dir, i), 'r') as f:
                    lines = f.readlines()
                for line in lines:
                    data = line.split()
                    if len(data) == 0:
                        continue
                    elif data[0] == 'Element':
                        ele_name[data[1]] = i
                        break
                    else:
                        continue
            for atoms_list_name in atoms_list:
                basis.append(join(basis_dir, ele_name[atoms_list_name]))

        else:
            raise ValueError("Please sure what you do!!!")
        
        return basis


def write_input_stru_core(stru=None,
                          directory="./",
                          filename="STRU",
                          potential=None,
                          pseudo_dir="./",
                          basis=None,
                          basis_dir="./",
                          coordinates_type="Cartesian",
                          atoms_list=None,
                          atoms_position=None, 
                          atoms_masses=None,
                          atoms_magnetism=None,
                          fix=1):
    if not judge_exist_stru(stru):
        return "No input structure!"

    elif (atoms_list is None):
        return "Please set right atoms list"
    elif(atoms_position is None):
        return "Please set right atoms position"
    elif(atoms_masses is None):
        return "Please set right atoms masses"
    elif(atoms_magnetism is None):
        return "Please set right atoms magnetism"
    else:
        with open(join(directory, filename), 'w') as f:
            f.write('ATOMIC_SPECIES\n')
            for i in range(len(atoms_list)):
                if(not exists(basename(potential[i]))):
                    shutil.copyfile(potential[i],
                                    directory+"/"+basename(potential[i]))
                temp1 = ' ' * (4-len(atoms_list[i]))
                temp2 = ' ' * (14-len(str(atoms_masses[i])))
                atomic_species = (atoms_list[i] + temp1
                                  + str(atoms_masses[i]) + temp2
                                  + basename(potential[i]))

                f.write(atomic_species)
                f.write('\n')

            f.write('\n')
            f.write('NUMERICAL_ORBITAL\n')
            for i in range(len(atoms_list)):
                if(not exists(basename(basis[i]))):
                    shutil.copyfile(basis[i],
                                        directory+"/"+basename(basis[i]))
                f.write(basename(basis[i]))
                f.write('\n')

            f.write('\n')
            f.write('LATTICE_CONSTANT\n')
            f.write('1.889726125 \n')
            f.write('\n')

            f.write('LATTICE_VECTORS\n')
            for i in range(3):
                for j in range(3):
                    temp3 = str("{:0<12f}".format(
                            stru.get_cell()[i][j])) + ' ' * 3
                    f.write(temp3)
                    f.write('   ')
                f.write('\n')
            f.write('\n')

            f.write('ATOMIC_POSITIONS\n')
            f.write(coordinates_type)
            f.write('\n')
            f.write('\n')
            for i in range(len(atoms_list)):
                f.write(atoms_list[i])
                f.write('\n')
                f.write(str("{:0<12f}".format(atoms_magnetism[i])))
                # update 20201230
                f.write('\n')
                f.write(str(len(atoms_position[i])))
                f.write('\n')

                for j in range(len(atoms_position[i])):
                    temp4 = str("{:0<12f}".format(
                            atoms_position[i][j][0])) + ' ' * 3
                    temp5 = str("{:0<12f}".format(
                            atoms_position[i][j][1])) + ' ' * 3
                    temp6 = str("{:0<12f}".format(
                            atoms_position[i][j][2])) + ' ' * 3
                    sym_pos = (temp4 + temp5 + temp6 +
                               (str(fix) + '   ') * 3)
                    f.write(sym_pos)
                    f.write('\n')
                # f.write('\n\n')

        pb_information = {}
        pb_information['pseudo_dir'] = pseudo_dir
        pb_information['basis_dir'] = basis_dir
        pb_information['potential_name'] = potential
        pb_information['basis_name'] = basis
        return pb_information


def write_input_stru(stru=None,
                     pseudo_dir='./',
                     potential_name=None,
                     basis_dir='./',
                     basis_name=None,
                     fix=1,
                     filename='STRU',
                     directory='./',
                     coordinates_type='Cartesian',
                     spin=1,
                     **kwargs):

    if not judge_exist_stru(stru):
        return "No input structure!"

    else:
        (atoms_list,
         atoms_masses,
         atoms_position,
         atoms_magnetism) = read_ase_stru(stru, coordinates_type)

        potential = set_potential(atoms_list,
                                  pseudo_dir,
                                  potential_name)
        basis = set_basis(atoms_list,
                          basis_dir,
                          basis_name)
        if(spin==2):
            for i in range(len(atoms_list)):
                atoms_magnetism[i] = 1.0
            

        pb_information = write_input_stru_core(stru,
                                               directory,
                                               filename,
                                               potential,
                                               pseudo_dir,
                                               basis,
                                               basis_dir,
                                               coordinates_type,
                                               atoms_list,
                                               atoms_position,
                                               atoms_masses,
                                               atoms_magnetism,
                                               fix)

        return pb_information


def read_stru(filename='STRU',
              directory='./',
              ase=True,
              **kwargs):
    # Read structure information from abacus structure file
    try:
        f = open(join(directory, filename), 'r')
    except Exception:
        return "Failed to open 'STRU', Please Check!"
    else:
        lines = f.readlines()
        f.close()

    # initialize reading information
    temp = []
    for line in lines:
        line = line.strip()
        line = line.replace('\n', ' ')
        line = line.replace('\t', ' ')
        line = line.replace('//', ' ')
        line = line.replace('#', ' ')

        if len(line) != 0:
            temp.append(line)

    atom_species = 0
    for i in range(len(temp)):
        if temp[i] == 'NUMERICAL_ORBITAL':
            atom_species = i - 1
            break

    atom_symbol = []
    atom_mass = []
    atom_potential = []
    atom_number = []
    atom_magnetism = []
    atom_positions = []
    atom_fix = []

    # get symbol, mass, potential
    for i in range(1, atom_species+1):
        atom_symbol.append(temp[i].split()[0])
        atom_mass.append(float(temp[i].split()[1]))
        atom_potential.append(temp[i].split()[2])
        atom_number.append(0)
        atom_magnetism.append(0)
        atom_positions.append([])
        atom_fix.append([])

    # get basis
    atom_basis = []
    for i in range(atom_species+2, (atom_species+1) * 2):
        atom_basis.append(temp[i].split()[0])

    # get lattice
    atom_lattice_scale = float(temp[(atom_species+1) * 2 + 1].split()[0])
    atom_lattice = np.array(
            [[float(temp[(atom_species+1) * 2 + 3 + i].split()[:3][j])
                for j in range(3)] for i in range(3)])

    # get coordinates type
    atom_coor = temp[(atom_species + 1) * 2 + 7].split()[0]

    # get position,  atoms number, magnetism, fix
    for i in range(atom_species):
        pos_start = (atom_species + 1) * 2 + 8 + 3 * i
        for j in range(i):
            pos_start += atom_number[j]
        atom_it = atom_symbol.index(temp[pos_start].split()[0])
        atom_magnetism[atom_it] = float(temp[pos_start + 1].split()[0])
        atom_number[atom_it] = int(temp[pos_start + 2].split()[0])

        atom_positions[atom_it] = np.array(
                [[float(temp[pos_start + 3 + i].split()[:3][j])
                    for j in range(3)] for i in range(atom_number[atom_it])])

        atom_fix[atom_it] = np.array(
                [[int(temp[pos_start + 3 + i].split()[3:6][j])
                    for j in range(3)]for i in range(atom_number[atom_it])])

    # Reset structure information and return results
    formula_symbol = ''
    formula_positions = []
    for i in range(atom_species):
        if atom_number[i] == 1:
            formula_symbol += atom_symbol[i]

        else:
            formula_symbol += atom_symbol[i] + str(atom_number[i])

        for j in range(atom_number[i]):
            formula_positions.append(atom_positions[i][j])

    formula_cell = atom_lattice * atom_lattice_scale * 0.529177210903

    if ase is True:
        if atom_coor == 'Direct':
            return Atoms(symbols=formula_symbol,
                         cell=formula_cell,
                         scaled_positions=formula_positions)

        elif atom_coor == 'Cartesian':
            return Atoms(symbols=formula_symbol,
                         cell=formula_cell,
                         positions=formula_positions)

        else:
            raise ValueError("atomic coordinate type is ERROR")

    else:
        return (formula_symbol,
                formula_cell,
                formula_positions,
                atom_potential,
                atom_basis)


if __name__ == "__main__":
    # print(potential_list())
    # print(basis_list())
    StruName = read_stru(filename='ABACUS_StruLi4Sn10 ',
                         directory='E:\Git\project\Structure',
                         ase=True)
    # print(StruName)
    print(write_input_stru(stru=StruName,
                           pseudo_dir='E:\Git\project\Potential\PotPotSG15',
                           potential_name=['Li_ONCV_PBE-1.0.upf',
                                           'Sn_ONCV_PBE-1.0.upf'],
                           basis_dir='E:\Git\project\Basis\BasSG15act',  
                           basis_name='SG15act',
                           fix=1,
                           filename='STRU',
                           directory='./',
                           coordinates_type='Cartesian'))
