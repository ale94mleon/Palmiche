#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from palmiche.utils import tools
from glob import glob
def wham(select_file_names = "coord[0-9]_selected.dat", OMP_NUM_THREADS = 12):
    select_file_names = sorted(glob(select_file_names))
    #wham_general = ["gmx wham -it tpr_files.dat -if pullf_files.dat -o all_coords_selected -hist hist_all_coords_selected -unit kJ"]
    wham_individual_coords = [f"export OMP_NUM_THREADS={OMP_NUM_THREADS}; gmx wham -ac -temp 303.15 -zprof0 4.5 -bins 300 -unit kJ -nBootstrap 200 -bs-method hist -is {s} -ix pullx_files.dat -it tpr_files.dat -o {s.split('.')[0]} -hist hist_{s.split('.')[0]} -oiact iact_{s.split('.')[0]}.xvg -bsres bsResult_{s.split('.')[0]}.xvg -bsprof bsProfs_{s.split('.')[0]}.xvg" for s in select_file_names]
    for cmd in wham_individual_coords: #wham_general + 
        tools.run(cmd)

if __name__ == '__main__':
    wham()