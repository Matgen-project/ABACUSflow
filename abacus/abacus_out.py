from __future__ import print_function
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 16:33:38 2018

Modified on Wed Jun 20 15:00:00 2018
@author: Shen Zhen-Xiong
"""

import subprocess
from os.path import join
import numpy as np
from ase.calculators.abacus.create_input import AbacusInput
from ase.calculators.calculator import Calculator, FileIOCalculator, all_changes  #Calculator


class Abacus(AbacusInput, FileIOCalculator):
    # Initialize parameters and get some information -START-
    name = 'Abacus'
    implemented_properties = ['energy', 'forces', 'fermi']

    default_parameters = dict(calculation='scf',
                              ecutwfc=50,
                              smearing='gaussian',
                              mixing_type='pulay-kerker',
                              basis_type='lcao',
                              gamma_only=1,
                              ks_solver="genelpa",
                              atom_file='STRU',
                              )

    def __init__(self,
                 restart=None,
                 ignore_bad_restart_file=False,
                 directory=None,
                 label='ase_rundir/ase_test',
                 atoms=None,
                 command=None,
                 log_file=None,
                 pseudo_dir='./',
                 potential_name=None,
                 basis_dir='./',
                 basis_name=None,
                 fix=1,
                 stru_filename='STRU',
                 coordinates_type="Cartesian",
                 **kwargs):

        self.species = None

        AbacusInput.__init__(self, restart)

        FileIOCalculator.__init__(self,
                                  restart,
                                  ignore_bad_restart_file,
                                  label,
                                  atoms,
                                  **kwargs)

        self.restart = restart
        self.pseudo_dir = pseudo_dir
        self.potential_name = potential_name
        self.basis_dir = basis_dir
        self.basis_name = basis_name
        self.fix = fix
        self.stru_filename = stru_filename
        self.coordinates_type = coordinates_type

        self.out_path = ''

        if directory is not None:
            self.directory = directory
        if log_file is not None:
            self.log_file = log_file

        AbacusInput.set(self, **self.parameters)
        AbacusInput.set(self, **kwargs)

    # Initialize parameters and get some information -END-

    def check_state(self, atoms):
        system_changes = FileIOCalculator.check_state(self, atoms)
        # Ignore boundary conditions:
        if 'pbc' in system_changes:
            system_changes.remove('pbc')
        return system_changes

    def initialize(self, atoms):
        numbers = atoms.get_atomic_numbers().copy()
        self.species = []
        for a, Z in enumerate(numbers):
            if Z not in self.species:
                self.species.append(Z)
        self.general_params["ntype"] = len(self.species)

    # Run abacus
    def calculate(self,
                  atoms=None,
                  properties=None,
                  system_changes=all_changes):
        FileIOCalculator.calculate(self,
                                   atoms,
                                   properties,
                                   system_changes)

    # Read results
    def read_results(self):
        a = AbacusInput()
        a.read_input_input(directory=self.directory)

        if self.pw_params['dr2'] is None:
            self.charge_density_error = float(1e-09)
        else:
            self.charge_density_error = float(self.pw_params['dr2'])

        if self.general_params['calculation'] in ['scf',
                                                  'relax',
                                                  'cell-relax',
                                                  'nscf',
                                                  'ienvelope',
                                                  'istate',
                                                  'test',
                                                  'md']:
            out_file = 'running_' + str(self.general_params['calculation']) + '.log'
        else:
            raise ValueError('Calculation parameters error')

        if self.general_params['suffix'] is None:
            self.out_path = join(self.directory, 'OUT.ABACUS/')
        else:

            self.out_path = join(self.directory, 'OUT.%s/' % str(self.general_params['suffix']))

        self.out_log_file = join(self.out_path, out_file)

        f = open(self.out_log_file, 'r')
        lines = f.readlines()
        f.close()

        n = 0
        number_atoms = 0
        force_number = 0
        force_last = []
        force_all = []
        fermi_energy = None
        final_total_energy = None

        for line in lines:
            if line.find('TOTAL ATOM NUMBER') != -1:
                number_atoms = int(line.split(' = ')[1])
            if line.find('final etot is') != -1:
                import re
                final_total_energy = re.findall(r'[-+]?\d+\.?\d*[eE]?[-+]?\d+', line)
                final_total_energy = float(final_total_energy[0])
            if line.find('EFERMI = ') != -1:
                fermi_energy = float(line.split()[2])
            if line.find('TOTAL-FORCE') != -1:
                for a in range(number_atoms):
                    force_all.append(
                              [float(data) for data in lines[n + 4 + a].split()[1:4]])
                force_number = force_number + 1
            n = n + 1

        force_all = np.array(force_all)
        force_last = force_all[-1 - number_atoms:]
        force_last = np.array(force_last)

        self.results['energy'] = final_total_energy
        self.results['fermi'] = fermi_energy
        self.results['forces'] = force_last

        return self.results

    def run(self):
        with open(self.log_file, 'a') as f:
            run = subprocess.Popen(self.command,
                                   stderr=f,
                                   stdin=f,
                                   stdout=f,
                                   cwd=self.directory,
                                   shell=True)
            return run.communicate()

    def get_fermi_level(self):
        return self.results['fermi']

    """
    def get_potential_energy(self, atoms):
        return self.get_property('energy', atoms)

    def get_forces(self, atoms):
        return self.get_property('forces', atoms)

    def get_property(self, name, atoms = None, allow_calculation = True):
        if atoms is None:
            atoms = self.atoms
            system_changes = []
        else:
            system_changes = self.check_state(atoms)
            if system_changes:
                self.reset()
        if name not in self.results:
            if not allow_calculation:
                return None
            self.calculate(atoms)
        result = self.results[name]
        return result
    """


if __name__ == "__main__":
    pass
