#!/bin/bash
echo 'start'
python modify_id.py sp &  
python modify_id.py dos_plotter &  
python modify_id.py bs_plotter &  
python modify_id.py formation_energy &  
python modify_id.py kpoints_data &  
python modify_id.py magnetization_property &  
python modify_id.py optimized_structure_cif &  
python modify_id.py pdf &  
python modify_id.py unoptimized_poscar &  
python modify_id.py xrd &  

wait
echo 'done!'
