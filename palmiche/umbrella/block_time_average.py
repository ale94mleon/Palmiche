#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from palmiche.utils import xvg, tools
import matplotlib.pyplot as plt
import os
import tqdm
import multiprocessing as mp
import copy
import re
def cmdwham(it = "tpr_files.dat",
            if_ = "pullf_files.dat",
            is_ = "coordsel.dat",
            unit = "kJ",
            temperature = 303.15,
            out = "out.xvg",
            hist = "hist.xvg",
            nice = None,
            b = None,
            e = None):
    string = f"gmx wham -it {it} -if {if_} -unit {unit} -is {is_} -temp {temperature} -o {out} -hist {hist}"
    if b: string+= f" -b {b}"
    if e: string+= f" -e {e}"
    if nice: string+= f" -nice {nice}"
    string += ''
    tools.run(string)
def cmdwham_star(keywords):
    return cmdwham(**keywords)

def pmf(time, NumBlocks, out_folder = 'BLOCKS', figsize = (16,9), figname = 'PMF_time_block_average.pdf', plot = True, wham = True, wham_jobs = 1, **keywords):
    f"""Give the time block average PMF

    Args:
        time (int): time in ps
        NumBlocks (int): number of blocks to analysis.
        out_folder (str, optional): Where aal the xvg will be saved. Defaults to 'BLOCKS'.
        figsize (tuple, optional): The size of the figure to matplotlib. Defaults to (10,9).
        figname (str, optional): The name of the figure to export. Defaults to PMF_time_block_average.pdf.
        plot (bool, optional): If True, a figure will be created with name {figname}. Defaults to False.
        wham (bool, optional): If True, gmx wham will be executed. Defaults to False.\
        wham_jobs (int, optional): How many wham jobs are launched in parallel. Defaults to 1.
        **keywords: All the necessary options of the function cmdwham.

    """
    intervals = range(0, time + int(time/NumBlocks), int(time/NumBlocks))
    tools.makedirs(out_folder)
    xvg_label_paths = []
    list_keywords = []
    for (i, init) in enumerate(intervals):
        try:
            keywords['b'] = init
            keywords['e'] = intervals[i+1]
            keywords['hist'] = os.path.join(out_folder, f'hist_{init}-{intervals[i+1]}.xvg')
            keywords['out'] = os.path.join(out_folder, f'pmf_{init}-{intervals[i+1]}.xvg')
            xvg_label_paths.append((f'{init}-{intervals[i+1]}', keywords['out']))
            list_keywords.append(copy.deepcopy(keywords))   
        except:
            pass
    if wham:
        if wham_jobs:
            jobs = wham_jobs
        else:
            jobs = mp.cpu_count()
        pool = mp.Pool(jobs)
        print(f"Generating PMF with gmx wham. It could take a long time. {jobs} job(s) was(were) launched (in parallel).")
        for i in tqdm.tqdm(pool.imap_unordered(cmdwham_star, list_keywords), total=len(list_keywords)):
            pass
        pool.close()
        print("Done!")

    whole_data = dict()
    for (label, path) in xvg_label_paths:
        whole_data[label] = xvg.XVG(path).data
    
    if plot:       
        fig, ax = plt.subplots(figsize = figsize)
        
        for key in whole_data:
            ax.plot(whole_data[key][:,0], whole_data[key][:,1], label=f"{key} ps")

        ax.legend(loc="lower right")
        ax.set(title = f"Time Block Average",
        xlabel = r"$\xi$",
        ylabel = "PMF [kJ/mol]")
        fig.savefig(figname)

    return whole_data

def multi_pmf(time, NumBlocks, path = '.', pattern = 'coord[0-5]_selected.dat', figsize = (16,9), figname = 'PMF_time_block_average.pdf', wham = True, wham_jobs = 1):
    linestyle_tuple = [
     ('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 1))),
     ('densely dotted',        (0, (1, 1))),

     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]
    
    pattern = re.compile(pattern)
    dats = [dat for dat in tools.list_if_file(path) if pattern.match(dat)]
    fig, ax = plt.subplots(figsize = figsize)
    if len(dats) > 1:
        fig_i, ax_i = plt.subplots(len(dats), 1,figsize = figsize)

    for (i,dat) in enumerate(dats):
        
        name = dat.split('.')[0]
        print(f'Calculating for {name}...')
        pmfData = pmf(time,
                        NumBlocks,
                        out_folder = f'BLOCKS_{name}',
                        plot = False,
                        wham_jobs = wham_jobs,
                        wham = wham,
                        is_ = dat)
        for (j, key) in enumerate(pmfData):
            try:
                linestyle = linestyle_tuple[j][1]
            except:
                linestyle = linestyle_tuple[j-len(pmfData)-1][1]
            ax.plot(pmfData[key][:,0], pmfData[key][:,1], label=f"{name}_{key} ps", linestyle = linestyle)
            if len(dats) >1:
                ax_i[i].plot(pmfData[key][:,0], pmfData[key][:,1], label=f"{name}_{key} ps")

        ax_i[i].legend(loc="lower right")
        ax_i[i].set(title = f"Time Block Average",
        xlabel = r"$\xi$",
        ylabel = "PMF [kJ/mol]",
        ylim = (-10,150))
        if i != len(dats) -1:
            ax_i[i].tick_params('x', labelbottom=False)



        


    ax.legend(loc="lower right")
    ax.set(title = f"Time Block Average",
    xlabel = r"$\xi$",
    ylabel = "PMF [kJ/mol]")

    fig.savefig(figname)
    fig_i.subplots_adjust(hspace=.5)
    fig_i.savefig(f"Split_{figname}", bbox_inches="tight") 

if __name__ == '__main__':
    pass

