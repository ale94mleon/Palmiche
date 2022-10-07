#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from palmiche.utils import tools, xvg
import argparse
import multiprocessing as mp
import os, re, tqdm
from scipy import integrate
import numpy as np


def gmx_analyze(f:str, ac:str = 'autocorr.xvg', overwrite:bool = False, reaction_coordinates:list = None, **kwargs):
    """Execute gmx analyze command

    Parameters
    ----------
    f : str
        The path to the xvg file with the format (time, value)
    ac : str, optional
        The name of the output for the autocorrelation calculations, by default 'autocorr.xvg'
    overwrite : bool, optional
        If True and if the ac file exists, this file will be overwrite,
        if False it will not perform any new calculation, by default False
    reaction_coordinates : list, optional
        This is a list for the reaction coordinate to take into account.
        If f has six columns, the first one will be considered time and the rest the reaction coordinates.
        If we would like to calculate for the second and third columns, then we must provided the list [2,3]
        by default None
    kwargs : Extra parameters will be passed directly to gmx analyze.
    """
    print(ac)
    if reaction_coordinates:
        columns = [0] + reaction_coordinates # we must pug the time
        f_dirname = os.path.dirname(f)
        f_name = os.path.basename(f).split('.')[0]
        f2work = os.path.join(f_dirname, f"{f_name}_acf_processed.xvg")
        f_xvg = xvg.XVG(f)
        f_xvg.data = f_xvg.data[:,columns]
        f_xvg.write(f2work)
    else:
        f2work = f

    if overwrite:
        cmd = 'export GMX_MAXBACKUP=-1; '
    else:
        cmd = ''
    cmd += f'gmx analyze -f {f2work} -ac {ac} '
    for key in kwargs:
        cmd += f'-{key} {kwargs[key]} '
    if not os.path.isfile(ac) or overwrite:
        tools.run(cmd)

def ___gmx_analyze_wrapper___(kwargs):
    """A simple wrapper to use parallelization of gmx_analyze.
    """
    gmx_analyze(**kwargs)

def IntegralOfACF(windows_path:str,  ac:str = 'autocorr.xvg', integrator:str = 'trapezoid', njobs:int = 1):
    """This will calculate the integral of the autocorrelation function

    Parameters
    ----------
    windows_path : str
        The path were the windows directories are with the format of of five integers. E.g. for window 1: 00001
    ac : str, optional
        The name of the output for the autocorrelation calculations, by default 'autocorr.xvg'
    integrator : str, optional
        Integration method. Are currently supported:
            trapezoid --> scipy.integrate.trapezoid
            trapz --> numpy.trapz (Default)
            simps --> scipy.integrate.simpson
    njobs : int, optional
        The number of jobs for the parallelization, by default 1

    Returns
    -------
    np.array of shape (number of windows X number of reaction coordinates)
    and a txt file with the saved dat in windows_path directory.

    Raises
    ------
    ValueError
        If an invalid integrator is provided.
    """
    pattern = re.compile("^\d{5}$")
    individual_paths = sorted([os.path.join(windows_path, window) for window in tools.list_if_dir(windows_path) if pattern.match(window)])

    # Selecting the integrator method
    if integrator == 'trapezoid':
        integrator = integrate.trapz
    elif integrator == 'trapz':
        integrator = np.trapz
    elif integrator == 'simpson':
        integrator = integrate.simps
    else:
        raise ValueError(f"The integration method '{integrator}' is not available.\n"\
            "Current possible selections are:\n"\
            "   -> trapezoid\n"\
            "   -> trapz\n"\
            "   -> simpson (Default)")

    # Getting all the data
    all_data = []
    for window in individual_paths:
        multidata = xvg.XVG(os.path.join(window, ac)).multidata

        data_on_window = multidata[0]
        for data in multidata[1:]:
            data_on_window = np.insert(data_on_window,-1, data[:,1], axis = -1)
        all_data.append(data_on_window)
    #all data shape = (number of windows, number of time points, number or reaction coordinates +1 [time is the first column])
    all_data = np.stack(all_data, axis = 0)

    # Get the IACFs per windows
    args_per_window = []
    for i in range(all_data.shape[0]):
        args_per_window.append((all_data[i], integrator))

    pool = mp.Pool(min(njobs, mp.cpu_count()))
    integrals_per_window = pool.starmap(___IntegralOfACF_chunk___, args_per_window)
    pool.close()
    integrals_per_window = np.asarray(integrals_per_window)
    # Saving data
    np.save(os.path.join(windows_path, 'IACFs_per_windows.npy'), integrals_per_window)

    # Get the integral on the time projection
    # Rearrange the data
    rearranged_data = []
    for j in range(all_data.shape[1]):
        rearranged_data.append(all_data[:,j,:])
    rearranged_data = np.asarray(rearranged_data)

    args_per_time = []
    for i in range(rearranged_data.shape[0]):
        args_per_time.append((rearranged_data[i], integrator))

    pool = mp.Pool(min(njobs, mp.cpu_count()))
    integrals_per_time = pool.starmap(___IntegralOfACF_chunk___, args_per_time)
    pool.close()
    integrals_per_time = np.asarray(integrals_per_time)
    # Saving data
    np.save(os.path.join(windows_path, 'IACFs_per_time.npy'), integrals_per_time)


    return integrals_per_window, integrals_per_time


def ___IntegralOfACF_chunk___(data, integrator):
    """Part of the code of IntegralOfACF
    """
    # index = int(os.path.basename(window))
    integrals = []
    for j in range(1, data.shape[1]):
        integrals.append(integrator(data[:,j], data[:,0]))
    return integrals

import numpy as np
import matplotlib.pyplot as plt

def custom_subplots(coords:list,x_size:int = 6):
    """This will create a grid plot where axes[3n+1,0] will be used for
    the contour plot, axes[3n,0] for the integral along the "X" and
    axes[3n+1,1] for the integral along the "Y"

    Parameters
    ----------
    coords : list
        A list of name of coordinates (could be numbers). The string Hello will be interpreted as 5 coordinates ['H', 'e', ...]
        The proper way is ['Hello']
    x_size : int, optional
        the size of the figure along "X", by default 6

    Returns
    -------
    tuples
        (fig, axes)
    """
    x_size = 6
    y_size = (len(coords))*x_size

    fig = plt.figure(figsize=(x_size, y_size)) # , layout="constrained"
    # Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
    # the size of the marginal axes and the main axes in both directions.
    # Also adjust the subplot parameters for a square plot.

    gs = fig.add_gridspec(3*len(coords), 2, width_ratios=(4, 1), height_ratios=len(coords)*[1, 4, 0.5], # 1, 4, 0.5
                        left=0.1, right=0.9, bottom=0.1, top=0.9,
                        wspace=0.05, hspace=0.05)
    # Create the Axes.
    axes = np.empty([3*len(coords), 2], dtype=object)
    for i in range(len(coords)):
        axes[3*i+1,0] = fig.add_subplot(gs[3*i+1,0])
        axes[3*i,0] = fig.add_subplot(gs[3*i,0])#, sharex=axes[2*i+1,0])
        axes[3*i,0].set_xticklabels([])
        axes[3*i+1,1] = fig.add_subplot(gs[3*i+1,1])#, sharey=axes[2*i+1,0])
        axes[3*i+1,1].set_yticklabels([])

        axes[3*i,0].set(
            title = f'Coord {coords[i]}'
        )
    return (fig, axes)

fig, axes = custom_subplots('abc')


def plot(windows_path, ac):
    pattern = re.compile("^\d{5}$")
    individual_paths = sorted([os.path.join(windows_path, window) for window in tools.list_if_dir(windows_path) if pattern.match(window)])

    all_data = []
    for window in individual_paths:
        multidata = xvg.XVG(os.path.join(window, ac)).multidata

        data_on_window = multidata[0]
        for data in multidata[1:]:
            data_on_window = np.insert(data_on_window,-1, data[:,1], axis = -1)
        all_data.append(data_on_window)
    #all data shape = (number of windows, number of time points, number or reaction coordinates +1 [time is the first column])
    all_data = np.stack(all_data, axis = 0)
    # Getting the levels
    minimum = -0.5 # all_data[:,:,1:].min() # -1
    maximum = 1 # all_data[:,:,1:].max() # 1
    levels = np.linspace(minimum, maximum, 10)



    cm = plt.get_cmap('viridis')#gist_rainbow viridis

    fig, axes = plt.subplots(nrows =  all_data.shape[2]-1, ncols = 1, figsize = (5,10), constrained_layout=True)
    plt.setp(axes.flat, ylabel="Window ID" , xlabel= "Time [ps]")



    # for ax, col in zip(axes[0], ['System']):
    #     ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
    #                 xycoords='axes fraction', textcoords='offset points',
    #                 size='large', ha='center', va='baseline', fontsize = 30)

    # for ax, row in zip(axes[:,0], range(all_data.shape[0])):
    #     ax.annotate(str(row), xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
    #                 xycoords=ax.yaxis.label, textcoords='offset points',
    #                 size='large', ha='right', va='center', fontsize = 30)

    X, Y = np.meshgrid(range(1, all_data.shape[0]+1),all_data[0,:,0])
    X, Y = np.meshgrid(all_data[0,:,0], range(1, all_data.shape[0]+1))

    for i in range(1, all_data.shape[-1]):
        axes[i-1].set(
            title = f'coord {i}'
        )
        Z = np.reshape(all_data[:,:,i], X.shape)
        CS = axes[i-1].contourf(X, Y, Z, cmap = cm, levels = levels)

    cbar = fig.colorbar(CS, ax=axes)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('ACF', rotation=270, fontsize = 15)


    fig.savefig("test.png", bbox_inches="tight")

def main(windows_path, pullx, ac = 'autocorr.xvg', overwrite = False, reaction_coordinates = None,  integrator:str = 'trapezoid', njobs = 1, **kwargs):

    pattern = re.compile("^\d{5}$")
    individual_paths = sorted([os.path.join(windows_path, window) for window in tools.list_if_dir(windows_path) if pattern.match(window)])

    ListOfKwargs = []
    extra_kwargs = kwargs.copy()
    for window in individual_paths:
        main_kwargs ={
            'f':os.path.join(window, pullx),
            'ac':os.path.join(window, ac),
            'overwrite':overwrite,
            'reaction_coordinates': reaction_coordinates,
            }
        if kwargs:
            extra_kwargs = kwargs.copy()
            extra_kwargs.update(main_kwargs)
            ListOfKwargs.append(extra_kwargs)
        else:
            ListOfKwargs.append(main_kwargs)

    print('Calculating ACFs..')
    pool = mp.Pool(min(njobs, mp.cpu_count()))
    pool.map(___gmx_analyze_wrapper___, ListOfKwargs)
    # pool.imap(gmx_analyze_wrapper, ListOfKwargs)
    pool.close()
    print('Done!')

    print('Calculating integral of ACFs per windows...')
    IntegralOfACF(windows_path = windows_path,  ac = ac, integrator =integrator, njobs=njobs)
    print('Done!')









if __name__ == '__main__':
    windows_path1 = '/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_O/7e27/sys_MMV007839_Cell_891_SP_param/windows'
    windows_path2 = '/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_Q/7e27/sys_MMV007839_Cell_891_SP_param/windows'
    windows_path3 = '/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_A/7e27/sys_MMV007839_Cell_891_SP_param/windows'

    ww = [windows_path3]
    # main(
    #     windows_path = windows_path,
    #     pullx='production_pullx.xvg',
    #     ac = 'autocorr.xvg',
    #     overwrite=True,
    #     reaction_coordinates=[1, 2, 3, 4, 5],
    #     njobs=12
    # )
    # my_dict = {'f': '/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_A/7e27/sys_MMV007839_Cell_891_SP_param/windows/00020/production_pullx.xvg', 'ac': '/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_A/7e27/sys_MMV007839_Cell_891_SP_param/windows/00020/autocorr.xvg', 'overwrite': False}
    # # gmx_analyze_wrapper(my_dict)
    # autocorr_xvg = xvg.XVG('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_A/7e27/sys_MMV007839_Cell_891_SP_param/windows/00020/autocorr.xvg')
    # print(autocorr_xvg.multidata[0])
    # print(len(autocorr_xvg.multidata))

    for i, w in enumerate(ww):
        a, b = IntegralOfACF(windows_path = w, njobs=12)
        print(a.shape, b.shape)
