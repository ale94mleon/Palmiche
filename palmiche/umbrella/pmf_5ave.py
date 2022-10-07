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
import re, os, tempfile, tqdm
from palmiche.tools import tools, xvg
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
# I will deffine my cero potential as the potential of the system when the ligand is in the bulk without any interaction (the last poitn in the PMF)

def cmdwham(it = "tpr_files.dat",
            if_ = "pullf_files.dat",
            is_ = "coordsel.dat",
            unit = "kJ",
            out = 'data.xvg'):
    if os.path.exists(is_):
        tools.run(f"export GMX_MAXBACKUP=-1;gmx wham -it {it} -if {if_} -unit {unit} -is {is_} -o {out} -nice 19")
    else:
        tools.run(f"export GMX_MAXBACKUP=-1;gmx wham -it {it} -if {if_} -unit {unit} -o {out}-nice 19")


def pmf_5ave(path, pattern = 'coord[1-5]_selected.dat', plot = False, wham = False, **keywords):
    if wham: tmpout = tempfile.NamedTemporaryFile(suffix='.xvg')
    cwd = os.getcwdb()
    path = os.path.abspath(path)
    pattern = re.compile(pattern)
    reaction_coord = []
    pmf = []
    os.chdir(path)
    print("Generating PMF:")
    for item in tqdm.tqdm([i for i in tools.list_if_file(path) if pattern.match(i)]):
        if wham:

            cmdwham(is_ = item, out = tmpout.name)
            data = xvg.XVG(tmpout.name).data
        else:
            data = xvg.XVG(item).data
        reaction_coord.append(data[:,0])
        pmf.append(data[:,1])
  

    data = np.asarray([
                    np.mean(reaction_coord, axis = 0),
                    stats.sem(reaction_coord, axis = 0, ddof=1),
                    np.mean(pmf, axis = 0),
                    stats.sem(pmf, axis = 0, ddof=1)
                    ]).transpose()
    #data = np.random.randint(0,30,size = data.shape)
    if plot:
        X = data[:,0]
        Xerror = data[:,1]
        Y = data[:,2]
        Yerror = data[:,3]
        fig, ax = plt.subplots(**keywords)

        ax.errorbar(X, Y,xerr = Xerror)
        ax.fill_between(X, Y - Yerror, Y + Yerror, alpha = 0.2)

        ax.set(title = f"Mean PMF using the five channels [Kj/mol]",
        xlabel = r"$\xi$",
        ylabel = "PMF [kJ/mol]")
        plt.show()
        fig.savefig("PMF_5ave.svg")
    os.chdir(cwd)
    return data

def main(start_with = 'umbrella', pattern = 'coord[1-5]_selected.dat', wham = False, **keywords):
    cwd = os.getcwd()
    paths = sorted([os.path.abspath(os.path.join(p, "7e27", "sys_MMV007839_Cell_891_SP_param", "windows")) for p in tools.list_if_dir(".") if p.startswith(start_with)])
    string = f"{'simulation':<100}{'Y_error_mean':<10}\n"
    fig, ax = plt.subplots(**keywords)
    
    for path in paths:
        print(path)
        data = pmf_5ave(path, pattern, plot = False, wham = wham)
        for s in path.split('/'):
            if s.startswith(start_with):
                label = s.split(start_with)[-1].strip()
        X = data[:,0]
        Xerror = data[:,1]
        Y = data[:,2]
        Yerror = data[:,3]

        ope = np.round(np.mean(Yerror),2)
        string += f"{label:<100}{ope:<10}\n"

        ax.plot(X, Y, label=label)
        ax.fill_between(X, Y - Yerror, Y + Yerror, alpha = 0.2)
        
    os.chdir(cwd)
    print(string)
    
    ax.legend(loc="lower right")
    ax.set(title = "Different PMF parameters. Using Average.",
    xlabel = r"$\xi$",
    ylabel = "PMF [kJ/mol]")
    fig.savefig("PMF_resume_5ave.svg")
    fig.show()
    
    with open("PMF_error_5ave.txt", "w") as f:
        f.write(string)
