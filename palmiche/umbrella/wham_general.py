#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : 2020-2023
Author        : Alejandro Martínez León
Mail          : [alejandro.martinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""
from palmiche.utils import tools
def wham():
    ligands = [0,1,2,3,4]
    number_windows = 87
    number_coords = 3*len(ligands)#THis case I use 3 because I use 2 pull restrain peer coord. orientation and loss of channel entrance 
    select_file_names = []
    for i in range(len(ligands) + 1):
        if i < len(ligands):
            row = tools.zerolistmaker(number_coords)
            row[i] = 1
            select_file_names.append(f"coord{i+1}_selected.dat")
            with open(f"coord{i+1}_selected.dat", "w") as select:    
                for j in range(number_windows):
                    select.write(f"{' '.join([str(r) for r in row])}\n")
        else:
            row = tools.zerolistmaker(number_coords)
            row[:len(ligands)] = len(ligands)*[1]
            select_file_names.append(f"coord0_selected.dat")
            with open(f"coord0_selected.dat", "w") as select:    
                for j in range(number_windows):
                    select.write(f"{' '.join([str(r) for r in row])}\n")

                               
         

    wham_individual_coords = [f"gmx wham -ac -temp 303.15 -zprof0 4.5 -bins 300 -unit kJ -nBootstrap 200 -bs-method hist -is {s} -ix pullx_files.dat -it tpr_files.dat -o {s.split('.')[0]} -hist hist_{s.split('.')[0]} -oiact iact_{s.split('.')[0]}.xvg -bsres bsResult_{s.split('.')[0]}.xvg -bsprof bsProfs_{s.split('.')[0]}.xvg" for s in select_file_names]

    for cmd in wham_individual_coords:#wham_general + wham_individual_coords:
        tools.run(cmd)
if __name__ == '__main__':
    wham()