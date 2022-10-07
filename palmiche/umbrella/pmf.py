#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from palmiche.utils import xvg, tools
import matplotlib.pyplot as plt
import os
import numpy as np

# I will deffine my cero potential as the potential of the system when the ligand is in the bulk without any interaction (the last poitn in the PMF)

def cmdwham(it = "tpr_files.dat",
            if_ = "pullf_files.dat",
            is_ = "coordsel.dat",
            unit = "kJ",
            nBootstrap = 100,
            bs_method = "b-hist",
            bsres = "bsResult.xvg",
            bsprof = "bsProfs.xvg"):
    if os.path.exists(is_):
        tools.run(f"gmx wham -it {it} -if {if_} -unit {unit} -is {is_} -nBootstrap {nBootstrap} -bs-method {bs_method} -bsres {bsres} -bsprof {bsprof} -nice 19")
    else:
        tools.run(f"gmx wham -it {it} -if {if_} -unit {unit} -nBootstrap {nBootstrap} -bs-method {bs_method} -bsres {bsres} -bsprof {bsprof} -nice 19")

def main(wham = True, nBootstrap = 100, figsize = (16,9), start_with = 'umbrella'):
    cwd = os.getcwd()
    paths = sorted([os.path.abspath(os.path.join(p, "7e27", "sys_MMV007839_Cell_891_SP_param", "windows")) for p in tools.list_if_dir(".") if p.startswith(start_with)])
    string = f"{'simulation':<100}{'Y_error_mean':<10}\n"
    fig, ax = plt.subplots(figsize = figsize)
    for path in paths:

        os.chdir(path)
        if wham:
            cmdwham(nBootstrap = nBootstrap, is_ = "coord0_selected.dat")
        data = xvg.XVG("bsResult.xvg").data
        for s in path.split('/'):
            if s.startswith(start_with):
                label = s.split(start_with)[-1].strip()
        X = data[:,0]
        Y = data[:,1] - data[-1:,1]
        Yerror = data[:,2]
        ope = np.round(np.mean(Yerror),2)
        string += f"{label:<100}{ope:<10}\n"

        ax.plot(X, Y, label=label)
        ax.fill_between(X, Y - Yerror, Y + Yerror, alpha = 0.2)
        
    os.chdir(cwd)
    print(string)
    
    ax.legend(loc="lower right")
    ax.set(title = f"Different PMF parameters. nBootstrap = {nBootstrap}",
    xlabel = r"$\xi$",
    ylabel = "PMF [kJ/mol]")
    fig.savefig("PMF_resume.svg")
    fig.show()
    
    with open("PMF_error_mean.txt", "w") as f:
        f.write(string)
