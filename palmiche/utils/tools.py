#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Fri Jun 25 21:00:10 2021
Author        : Alejandro Martínez León
Mail          : [alejandro.martinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

import os, subprocess, shutil, datetime, tempfile, time, tqdm
from collections import OrderedDict
import numpy as np
import pandas as pd
import multiprocessing as mp
from glob import glob
from kmodes.kmodes import KModes

import prolif as plf
import MDAnalysis as mda
from MDAnalysis.topology.guessers import guess_types

from rdkit.Chem import DataStructs, AllChem
from rdkit.DataStructs import ExplicitBitVect
from rdkit import Chem


from sklearn.cluster import KMeans
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from scipy.cluster.vq import vq
from scipy.interpolate import interp1d
from palmiche.utils import xvg
from palmiche.home import home
from typing import Union
import tarfile

PathLike = Union[os.PathLike, str, bytes]


#=======================================================================================

#                          Miscellanea tools

#=======================================================================================
class CTE:
    """Some importnat physical constants
    """
    Kb = 8.314462618E-3 #kJ/(mol⋅K) (kNA)
    R = 8.314462618E-3 #kJ/(mol⋅K)

def get_palmiche_data(file:PathLike, out_dir:PathLike = '.') -> None:
    """Get data from the data directory of Palmiche


    Parameters
    ----------
    file : PathLike
        this is the path relative to the data directory of Palmiche. For example,
        To get the amber99sb-star-ildn force field, you must provided
        GROMACS.ff/amber99sb-star-ildn.ff.tar.gz
    out_dir : PathLike, optional
        Where the file will be decompress, by default '.'
    """
    fname = home(file)
    tar = tarfile.open(fname, "r:gz")
    tar.extractall(out_dir)
    tar.close()

def RSS(data, AlignPoint, ref_index_data = 0, NumbPoints = None, InterpolationKind = 'cubic'):
    """This function will return a list of Residues Squared Sums.

    Args:
        data (list): A list of np.arrays build in such a way that if reference = data[0]. Then for the reference: reference[:,0] give the independent variable and reference[:,1] the dependent one. NAd the same for the other entrances The first entrance of this list must be the reference
        AlignPoint (_type_): _description_
        NumbPoints (int, optional): _description_. Defaults to 1000.
        InterpolationKind (str, optional): _description_. Defaults to 'cubic'.

    Raises:
        ValueError: _description_

    Returns:
        _type_: _description_
    """
    # Example:
    # In [1]: from palmiche.utils import tools, xvg
    # In [2]: d0 = xvg.XVG('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_N_ST/7e27/sys_MMV007839_Cell_891_SP_param/windows/coord0_selected.xvg').data
    # In [3]: d1 = xvg.XVG('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_N_ST/7e27/sys_MMV007839_Cell_891_SP_param/windows/coord1_selected.xvg').data
    # In [4]: d2 = xvg.XVG('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_N_ST/7e27/sys_MMV007839_Cell_891_SP_param/windows/coord2_selected.xvg').data
    # In [5]: d3 = xvg.XVG('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_N_ST/7e27/sys_MMV007839_Cell_891_SP_param/windows/coord3_selected.xvg').data
    # In [6]: d4 = xvg.XVG('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_N_ST/7e27/sys_MMV007839_Cell_891_SP_param/windows/coord4_selected.xvg').data
    # In [7]: d5 = xvg.XVG('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_N_ST/7e27/sys_MMV007839_Cell_891_SP_param/windows/coord5_selected.xvg').data
    # In [8]: data = [d0, d1, d2, d3, d4, d5]
    # In [9]: AlignPoint = 4.5
    # In [10]: NumbPoints = 300
    # In [11]: print(tools.RSS(data,AlignPoint, NumbPoints = NumbPoints))
    # [0.0, 962.7853017103384, 4159.864180622066, 9456.312066256418, 2076.1281115807665, 4421.657283337928]

    # GEtting the number of points for the interpolation
    if not NumbPoints:
        NumbPoints = data[0].shape[0]
    # Setting the reference
    ref = data[ref_index_data]
    # Getting the interpolation data
    f_interp_data = [interp1d(ref[:,0], item[:,1], kind=InterpolationKind) for item in data]

    minimum = max([item[:,0].min() for item in data])
    maximum = min([item[:,0].max() for item in data])
    # In this way we ensure only take those values for which we already have points in both data set

    if  not minimum <=AlignPoint <= maximum:
        raise ValueError(f"AlignPoint = {AlignPoint} is not in the range of the data, which is : {(minimum, maximum)}")

    # Interpolating data nad having the same points (newx)
    newx = np.linspace(minimum, maximum, num=NumbPoints, endpoint=True)
    interp_data = [f_interp(newx) for f_interp  in f_interp_data]


    #Aligning all the data
    for i in range(len(data)):
        if i == ref_index_data:
            continue
        else:
            offset =  f_interp_data[i](AlignPoint) - f_interp_data[ref_index_data](AlignPoint)
            interp_data[i] -= offset
    # import matplotlib.pyplot as plt
    # for i in range(len(data)):
    #     plt.scatter(data[i][:,0], data[i][:,1], label = f'row_data{i}')
    #     plt.plot(newx, interp_data[i],'->', label = f'interp_data{i}', ms = 4)
    # plt.legend()
    # plt.show()
    # Getting RSS values
    RSS_values = []
    for i in range(len(interp_data)):
        RSS_values.append(((interp_data[i] - interp_data[ref_index_data])**2).sum())

    return RSS_values

def RUS(data, ref_index_data = 0, NumbPoints = None, InterpolationKind = 'cubic'):
    """
    RUS stands for Residues Unsigned Sums
    RUS = sum_{i}^{n} (M_i - R_i)
    Where M_i is the error value of the magnitude on i of the model and R_i for the reference.
    RUS < 0, The model presents lower errors than the reference (in general)
    RUS = 0, The model has the same errors as the reference or they cancell out
    RUS > 0, The model presents worst errors than the reference
    """
    # GEtting the number of points for the interpolation
    if not NumbPoints:
        NumbPoints = data[0].shape[0]
    # Setting the reference
    ref = data[ref_index_data]
    # Getting the interpolation data
    f_interp_data = [interp1d(ref[:,0], item[:,1], kind=InterpolationKind) for item in data]

    minimum = max([item[:,0].min() for item in data])
    maximum = min([item[:,0].max() for item in data])
    # In this way we ensure only take those values for which we already have points in both data set


    # Interpolating data nad having the same points (newx)
    newx = np.linspace(minimum, maximum, num=NumbPoints, endpoint=True)
    interp_data = [f_interp(newx) for f_interp  in f_interp_data]

    # import matplotlib.pyplot as plt
    # for i in range(len(data)):
    #     plt.scatter(data[i][:,0], data[i][:,1], label = f'row_data{i}')
    #     plt.plot(newx, interp_data[i],'->', label = f'interp_data{i}', ms = 4)
    # plt.legend()
    # plt.show()
    # Getting RSS values
    RUS_values = []
    for i in range(len(interp_data)):
        RUS_values.append(((interp_data[i] - interp_data[ref_index_data])).sum())

    return RUS_values

def sign_change_index(values: list) -> int:
    """This function check wheter the sign of the values change.
    For example, for the list [-3,-2,0,1,2], the return value will be 2.
    On value "0" the sig change.

    Parameters
    ----------
    values : list
        Iterable of numerical values

    Returns
    -------
    int
        The index. In case that the sign does not change or values has less than 2 elements,
        it will return None.
    """
    if len(values) <= 1:
        return None
    else:
        for i in range(1, len(values)):
            if values[i] * values[i-1] < 0:
                return i
    return None

def zerolistmaker(n):
    return [0]*n
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

    # >>> angle_between((1, 0, 0), (0, 1, 0))
    # 1.5707963267948966
    # >>> angle_between((1, 0, 0), (1, 0, 0))
    # 0.0
    # >>> angle_between((1, 0, 0), (-1, 0, 0))
    # 3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
def aovec(vectors, round = None):
    """Return the average oriented vector
    Calculate all the angles respect to the cartessina axis and average the angles,andthen return the unitary vector
    ang(vec, x) = mean(ang(vectors_i, x))
    ang(vec, y) = mean(ang(vectors_i, y))
    ang(vec, z) = mean(ang(vectors_i, z))
    vec = (cos(ang(vec, x); cos(ang(vec, y); cos(ang(vec, z))
    Args:
        vectors ([type]): [description]
        round ([type]): [description]
    Returns:
        [type]: [description]
    """
    vectors = [np.asarray(v) for v in vectors]
    #Dimension of the vectors, how many components
    vec = np.empty(vectors[0].shape)
    for i in range(vectors[0].shape[0]):
        zero_vector = np.zeros(vectors[0].shape)
        zero_vector[i] = 1
        vec[i] = np.cos(np.mean([angle_between(vector, zero_vector) for vector in vectors]))

    vec = unit_vector(vec)
    if round: vec = np.round(vec, round)
    return vec

def get_vec_COM(conf, tpr, ndx, group1, group2, round = None):
    """Take two groups of a index file and return the vector formed by the coordinates of group2, group1 in this order.py

    Args:
        conf (path): the configuration file [gro, pdb, xtc, etc..]
        tpr (path): tpr file, or gro, pdb, etc...
        ndx (path): index file, need to have the defined groups
        group1 (str): name of the first group in the index file.
        group2 (str): name of the second group in the index file.
        round (int or None): the specification to round.

    Returns:
        np.array: the unitary vector formed by group2, group1
    """
    with tempfile.NamedTemporaryFile(suffix=".xvg") as xvgtmp1:
        run(f"export GMX_MAXBACKUP=-1;echo {group1} | gmx traj -f {conf} -s {tpr} -n {ndx} -com -ox {xvgtmp1.name}")
        coord1 = xvg.XVG(xvgtmp1.name).data[0,1:]# First row from column 1 to the end, the column 0 is time.

    with tempfile.NamedTemporaryFile(suffix=".xvg") as xvgtmp2:
        run(f"export GMX_MAXBACKUP=-1;echo {group2} | gmx traj -f {conf} -s {tpr} -n {ndx} -com -ox {xvgtmp2.name}")
        coord2 = xvg.XVG(xvgtmp2.name).data[0,1:]# First row from column 1 to the end, the column 0 is time.
    vector = unit_vector(coord2 - coord1)
    if round: vector = np.round(vector, round)
    return vector


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print('%r  %2.2f ms' % \
                  (method.__name__, (te - ts) * 1000))
        return result
    return timed

#=======================================================================================

#                          Tools for execution

#=======================================================================================
def multi_run(commands, nPar, shell = True, executable = '/bin/bash'):
    """This will run as many runs as nPar.

    Args:
        commands (list): A list of string to be run in the specified shell.
        nPar (int): How many processes are running simultaneously
        shell (bool, optional): belongs to run(). Defaults to True.
        executable (str, optional): belongs to tun(). Defaults to '/bin/bash'.
        Popen (bool, optional): belongs to run(). Defaults to False.
    """

    tmp_cmds = []
    for cmd in commands:
        tmp_cmds.append(run(cmd, shell = shell, executable = executable, Popen = True))
        if len(tmp_cmds) >= nPar:
            print(f"{len(tmp_cmds)} running, now waiting...")
            for tmp_cmd in tmp_cmds:
                tmp_cmd.wait()
            tmp_cmds = []


def run(command, shell = True, executable = '/bin/bash', Popen = False):
    #Here I could make some modification in order that detect the operator system
    #NAd make the command compatible with the operator system
    #the function eval could be an option if some modification to the variable command
    #need to be done.... SOme flight ideas...

    if Popen:
        #In this case you could access the pid as: run.pid
        process = subprocess.Popen(command, shell = shell, executable = executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
    else:
        process = subprocess.run(command, shell = shell, executable = executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
        returncode = process.returncode
        if returncode != 0:
            print(f'Command {command} returned non-zero exit status {returncode}')
            raise RuntimeError(process.stderr)
    return process

class CHECK:
    """
    This class test the normal ending of a GROMACS simulation using the log file.
    E.g.
    let say that we have the md.log of a simulation
    # CHECK = CHECK('md.log')
    # if CHECK.performance:
    #     #The simulation had a good ending
    # else:
    #     #The simulations crashed.
    """

    def __init__(self, logfile, lisfile = None):
        self.logfile = logfile
        self.lisfile = lisfile
        self.performance = None  #in ns/day, this attribute could be used to check if the simulation ended
        self.__check_normal_end__()

    def __check_normal_end__(self):
        with open(self.logfile, 'r') as f:
            lines = f.readlines()
        for i in range(len(lines)):
            if 'Finished mdrun on rank' in lines[i]:
                try:
                    self.performance = float(lines[i - 1].split()[1])
                except:
                    pass
                break

    def estimated_end_datetime(self):
        """Read an md.lis looking for the last instance of the string with the estimation of the ended time.
                "step 24940000, will finish Tue Oct  5 10:38:53 2021"

        Args:
            self.lisfile (path): md.lis, the output of the gromacs execution.

        Returns:
            datetime: The estimated end of the simulation.
        """
        with open(self.lisfile, "r") as f:
            lines = f.readlines()
        datetime_string = None
        for line in reversed(lines):
            if "will finish" in line:
                datetime_string = line.split("will finish")[-1].strip()
                break

        if datetime_string:
            return datetime.datetime.strptime(datetime_string, "%a %b %d %H:%M:%S %Y")
        else:
            return None

    def get_startdatetime_from_log(self):
        with open(self.logfile, "r") as f:
            lines = f.readlines()
        datetime_string = None
        for line in reversed(lines):
            if 'Started mdrun on rank' in line: # Will take the last
                datetime_string = " ".join(line.split()[-5:])
                break

        if datetime_string:
            return datetime.datetime.strptime(datetime_string, "%a %b %d %H:%M:%S %Y")
        else:
            return None

    def daysdiff(self, startdatetime = None):
        """Take the difference between the estimated ended time of lisfile and datetime_start

        Args:
            lisfile (path): [description]
            datetime_start (datetime): A date time.

        Returns:
            float: Days' difference.
        """
        if not startdatetime:
            startdatetime = self.get_startdatetime_from_log()
            print("The starting data time of the log file was used becasue user did not provided a start_datetime")

        end_datetime = self.estimated_end_datetime()
        if startdatetime and end_datetime:
            return (end_datetime - startdatetime).total_seconds() / 86400
        else:
            return None

def checkrun(user = '$USER', partition = 'deflt'):
    """Check for the running process on the cluster.

    Returns:
        list of integers: The integers that identify the running process on the cluster.
    """
    stdout = run(f"squeue -p {partition} -u {user} --format=%.i").stdout.split()[1:]#The first line is the string "JOBID"

    return [int(item) for item in stdout]

def job_launch_list(job_path_list, shell="sbatch"):
    """
    Same as job_launch, but the path to the jobs are provided

    Args:
    job_path_list (_type_): _description_
    shell (str, optional): _description_. Defaults to "sbatch".

    Returns:
    _type_: _description_
    """
    cwd = os.getcwd()
    JOBIDs = []
    for job_path in job_path_list:
        path = os.path.dirname(job_path)
        script_name = os.path.basename(job_path)
        os.chdir(path)
        if shell == "sbatch":
            JOBIDs.append(int(run(f"{shell} --parsable {script_name}").stdout))
        else:
            try:
                JOBIDs.append(int(run(f"{shell} {script_name}").stdout))
            except:
                pass
    os.chdir(cwd)
    return JOBIDs

def launch_wait_check_repeat(partition, jobpaths, logfile, lisfile, jobsh_name, logfile_start_datetime = None):
    maximum_days = {
        'deflt': 2,
        'long': 5,
        'debug': 0.02,
        'medium': 2,
        'fat+': 2,
        'fat': 2,
        'int': 2,
        'gpu': 2,
        'gpu:test': 0.02,
        'gpu-hub':2,
        'uds-hub':2,
            }
    # partition = 'deflt'
    # jobs2launch = [full_path_jobsh]
    # logfile = "production.log"
    # lisfile = 'production.lis'
    # jobsh_name = 'job.sh'
    # First launch (normal)
    # I will check for the cpt file based on the logfile name
    cpt_file = os.path.splitext(logfile)[0] + '.cpt'
    JOBIDs = job_launch_list([os.path.join(path, jobsh_name) for path in jobpaths])
    while len(set(JOBIDs).difference(checkrun(partition=partition))) != len(JOBIDs):
        time.sleep(np.random.randint(20,60))
    # Relaunch if needed.
    jobs2relaunch = {}
    for path in jobpaths:
        check = CHECK(os.path.join(path,logfile), os.path.join(path,lisfile))
        # Check from the correct end
        print(f"The simulation {os.path.join(path,logfile)} presented a performance of {check.performance}")
        if not check.performance:
            if os.path.exists(cpt_file):
                if logfile_start_datetime:
                    tmp_check = CHECK(os.path.join(path,logfile_start_datetime))
                    daysdiff = check.daysdiff(tmp_check.get_startdatetime_from_log())
                else:
                    daysdiff = check.daysdiff()
                # Check if the estimated time is larger than the reservation time in the corresponded partition.
                if daysdiff:
                    if partition not in maximum_days:
                        maximum_days[partition] = 2
                    if daysdiff > maximum_days[partition]:
                        # Get how many time the simulation must be relaunch
                        repeat = round((daysdiff / maximum_days[partition]) - 1)
                        if repeat == 0:repeat = 1
                        # Add path as key and number of repetition as value
                        jobs2relaunch[path] = repeat
                    else:
                        print(f'We were not able to guess what is the problem with the run of the job {path}')
                else:
                    print(f'We were not able to guess what is the problem with the run of the job {path}')
            else:
                print(f'Unable to continue job {path}, {cpt_file} doe not exist or is unaccesible')
    # Only if there are some jobs to relaunch
    if jobs2relaunch:
        while sum(jobs2relaunch.values()):
            # Launch the jobs if repeat is greater than zero
            to_launch = [os.path.join(path, jobsh_name) for path in jobs2relaunch if jobs2relaunch[path] > 0]
            JOBIDs = job_launch_list(to_launch)
            while len(list(set(JOBIDs).difference(checkrun(partition=partition)))) != len(JOBIDs):
                time.sleep(np.random.randint(20,60))
            # Update the repetitions
            for path in jobs2relaunch:
                check = CHECK(os.path.join(path,logfile), os.path.join(path,lisfile))
                if check.performance:
                    jobs2relaunch[path] == 0
                if jobs2relaunch[path] > 0:
                    jobs2relaunch[path] -= 1
        print(f'Performance of the relaunched jobs:')
        for path in jobs2relaunch:
            check = CHECK(os.path.join(path,logfile), os.path.join(path,lisfile))
            print(f"The simulation {os.path.join(path,logfile)} presented a performance of {check.performance}")


def job_launch(shell="sbatch", script_name = "job.sh"):
    """


    Parameters
    ----------
    shell : TYPE, optional
        DESCRIPTION. The default is "sbatch".
    script_name : TYPE, optional
        DESCRIPTION. The default is "job.sh".
        If a regular expresion is provided
        then the function will execute the first ordered alphabetically. E.g:
        job.* was provided and there job.sh and job.bash. Then it will use job.bash.


    Returns
    -------
    JOBIDs : list of integers.
        DESCRIPTION. The ID of the launch in case of sbatch was used as shell

    """
    cwd = os.getcwd()
    JOBIDs = []
    for (root, _, _) in list(os.walk(cwd)):
        os.chdir(root)
        if glob(script_name):
            script_name = sorted(glob(script_name))[0]
            if shell == "sbatch":
                JOBIDs.append(int(run(f"{shell} --parsable {script_name}").stdout))
            else:
                try:
                    JOBIDs.append(int(run(f"{shell} {script_name}").stdout))
                except:
                    pass
    os.chdir(cwd)
    return JOBIDs


def backoff(file_path):
    """


    Parameters
    ----------
    file : TYPE string
        DESCRIPTION: The name or the path f     or the specific file

    Returns
    -------
    None.
    If the file already exist. it will made a back up to ./#{file}.{str(i)}#,
    Where i is an integer.
    """
    basname = os.path.basename(file_path)
    dirname = os.path.dirname(file_path)
    if os.path.exists(file_path):
        new_basname = basname
        i = 1
        while(os.path.exists(os.path.join(dirname, new_basname))):
            new_basname = f"./#{basname}.{str(i)}#"
            i += 1
        print(f"Back Off! I just backed up {file_path} to {os.path.join(dirname, new_basname)}")
        shutil.copy2(file_path, os.path.join(dirname, new_basname))


#=======================================================================================

#                          Tools for working with files

#=======================================================================================
def makedirs(path):
    if os.path.exists(path):
        pass
    else:
        os.makedirs(path,exist_ok=True)

def rm(pattern, r = False):
    """


    Parameters
    ----------
    patterns : string
        input-like Unix rm
    r : TYPE, bool
        DESCRIPTION. The default is False.
        If True delete also directories
    Returns
    -------
    None.

    """
    if glob(pattern):
        for p in glob(pattern):

            if r:
                if os.path.isdir(p):
                    shutil.rmtree(p)
                elif os.path.isfile(p):
                    os.remove(p)
            else:
                if os.path.isfile(p):
                    os.remove(p)
                elif os.path.isdir(p):
                    print(f"The directory '{p}' was not deleted, set the flag r to True")
    else:

        print(f"rm: '{pattern}' doesn't exist")


def cp(src, dest, r = False):
    """
    This function makes use of the possible multiple CPU of the machine.

    Parameters
    ----------
    src : TYPE: string
        DESCRIPTION:
        Source Path to a directory or a file. regular expresion are accepted. The librery glob is used for that.
    dest : TYPE: string
        DESCRIPTION:
        Destination Path
    r : TYPE, optional
        DESCRIPTION. The default is False:
        If True, Also directories will be copy. If not, and a directory was
        given as src, a Raise Exception will be printed
        Another Raise Exception will be printed if the destination path doesn't
        exist.

    Returns
    -------
    None.

    """
    src = os.path.abspath(src)
    dest = os.path.abspath(dest)
    #This scheme doesn't consider the possibilities of other item except files, dir or regualr expresions
    if glob(src) and os.path.exists(os.path.dirname(dest)):

        if os.path.isdir(src):
            if r == True:

                basename_src = os.path.basename(src)
                for root, dirs, files in os.walk(src):
                    path2copy = root.split(src)[1]
                    try:
                        if path2copy[0] == "/":
                            path2copy = path2copy[1:]
                    except: pass
                    path_dest = os.path.join(dest,basename_src,path2copy)

                    if dirs:
                        pool = mp.Pool(mp.cpu_count())
                        pool.map(makedirs, [os.path.join(path_dest, d) for d in dirs])
                        pool.close()
                    if files:
                        makedirs(path_dest)
                        pool = mp.Pool(mp.cpu_count())

                        pool.starmap(shutil.copy2, [(os.path.join(root, f),path_dest) for f in files])
                        pool.close()
            else:
                print(f"If you need to copy directories, set the 'r' flag to True. {src} is a directory")

        elif os.path.isfile(src):
            shutil.copy2(src,dest)
            #if not os.path.exists(dest):
             #   print(f"The file '{src} was not copy to '{dest}' becasue doesn't exist or is not accesible")


        else:#Here we are in regular expresions
            to_copy = [file for file in glob(src) if os.path.isfile(file)]
            pool = mp.Pool(mp.cpu_count())
            pool.starmap(shutil.copy2, [(pattern,dest) for pattern in to_copy])
            pool.close()

    elif not glob(src) and not os.path.exists(os.path.dirname(dest)):
        print(f"cp: neither {src}, nor {dest} exist.")
    elif not glob(src):
        print(f"cp: The source file {src} doesn't exist")
    else:
        print(f"cp: cannot create regular file '{dest}': Not a directory")



def mv(src, dest, r = False):
    cp(src, dest, r = r)
    rm(src,r=r)

def touch(filepath):
    open(filepath, "w").close()

def list_if_dir(path = '.'):
    return [item for item in os.listdir(path) if os.path.isdir(os.path.join(path, item))]

def list_if_file(path = '.'):
    return [item for item in os.listdir(path) if os.path.isfile(os.path.join(path, item))]

#=======================================================================================

#                      Tools for working with index and topology files

#=======================================================================================
def replace_conformer(mol: Chem.rdchem.Mol, ref_mol: Chem.rdchem.Mol, inplace: bool = True):
    """Will replace the conformation state of mol
    by the conformational state of ref_mol.
    Preserving the atom order of mol.

    Parameters
    ----------
    mol : Chem.rdchem.Mol
        The molecule to add the conformation
    ref_mol : Chem.rdchem.Mol
        The molecule with the reference conformation
    inplace : bool, optional
        If True mol will be modify inplace, if False a new instance will be return, by default True

    Returns
    -------
    Chem.rdchem.Mol
        The new instance of mol only if inplace = True, if not None will be return

    Raises
    ------
    ValueError
        If ref_mol does not have a valid conformational state.
    """
    if not inplace:
        from copy import deepcopy
        mol = deepcopy(mol)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(ref_mol), f"\n{Chem.MolToSmiles(mol)}\n{Chem.MolToSmiles(ref_mol)}"

    try:
         mol_conf = mol.GetConformer()
    except ValueError:
        AllChem.EmbedMolecule(mol)
        mol_conf = mol.GetConformer()

    try:
        ref_mol.GetConformer()
    except ValueError:
        raise ValueError('ref_mol does not have a conformer.')

    mcs_atoms_idxs = mol.GetSubstructMatch(ref_mol)
    for ref_mol_atom_idx, mol_atom_idx in enumerate(mcs_atoms_idxs):
        mol_conf.SetAtomPosition(mol_atom_idx, ref_mol.GetConformer().GetAtomPosition(ref_mol_atom_idx))
    if inplace:
        return None
    else:
        return mol

class Mol2:
    """Read a mol2 molecule and storge the corresponded RDKit molecule in
    the attribute 'mol'

    Example
    -------
    .. ipython:: python

        from palmiche.utils.tools import Mol2
        from palmiche import home
        mol2_path = home(os.path.join("sample", "test.mol2"))
        molecule = Mol2(mol2_path)
        for atom in molecule.mol.GetAtoms():
            # Two etra properties added to the atom class
            print(atom.GetProp('atom_type'), atom.GetDoubleProp('charge'))
    """
    def __init__(self, mol2_path:PathLike) -> None:
        self.mol2_path = mol2_path
        self.mol = None
        self._name_changed_ = False
        self.read()

    def read(self) -> None:
        """This function read the file and create
        the RDKit molecule and stores it in the attribute
        `mol`.

        """
        self.mol = Chem.MolFromMol2File(self.mol2_path, removeHs = False)
        if not self.mol:
            print('Trying with atom types change...')
            with open(self.mol2_path, 'r') as f:
                lines = f.readlines()
                atom_types = []
                user_formal_charges = []
                for i in range(len(lines)):
                    if lines[i].startswith('@<TRIPOS>ATOM'):
                        for j in range(i+1, len(lines)):
                            if lines[j].startswith('@'):
                                break
                            line_split = lines[j].split()
                            atom_types.append(line_split[5])
                            user_formal_charges.append(float(line_split[-1]))
                            line_split[5] = line_split[1].translate({ord(i): None for i in '0123456789'})
                            lines[j] = '\t'.join(line_split) + '\n'
                        break
            string = "".join(lines)
            self.mol = Chem.MolFromMol2Block(string, removeHs = False)
            if self.mol:
                for atom, atom_type, user_formal_charge in zip(self.mol.GetAtoms(), atom_types, user_formal_charges):
                    atom.SetDoubleProp('charge', float(user_formal_charge))
                    atom.SetProp('atom_type',atom_type)
                print("The mol was created succesfully!")
            else:
                print("\nIt is not possible to convert the mol2 to RDKit molecule.\n")

    def AddConformer(self, ref_mol:Chem.rdchem.Mol) -> None:
        """Add to self.mol the conformation of ref_mol.
        It is not needed that the molecules have the same atom
        indexes; `palmiche.utils.tools.replace_conformer` will
        be internally called.

        Parameters
        ----------
        ref_mol : Chem.rdchem.Mol
            The reference molecule with a conformational state.
        """
        replace_conformer(self.mol, ref_mol, inplace=True)

    def ChangeName(self, mol_name:str) -> None:
        """Change name of the molecule

        Parameters
        ----------
        mol_name : str
            Code name to the molecule. It is wise to use a code
            no larger than 3 characters.
        """
        self._name_changed_ = True
        self.mol.SetProp('_Name', mol_name)

    def write(self, out_file:PathLike) -> None:
        """Write the molecue in mol2 format.

        Parameters
        ----------
        out_file : PathLike
            The path to the file.
        """
        with open(self.mol2_path, 'r') as f:
            lines = f.readlines()

        for i in range(len(lines)):
            if lines[i].startswith('@<TRIPOS>ATOM'):
                atom_idx = 0
                for j in range(i+1, len(lines)):
                    if lines[j].startswith('@'):
                        break
                    line_split = lines[j].split()
                    line_split[2:5] = [coord for coord in self.mol.GetConformer().GetAtomPosition(atom_idx)]
                    if self._name_changed_:
                        line_split[7] = self.mol.GetProp('_Name')
                    for column, item in enumerate(line_split):
                        if column in [0,1,5,6,7]:
                            line_split[column] = f"{line_split[column]:>5}"
                        elif column in [2,3,4]:
                            line_split[2:5] = [float(coord) for coord in line_split[2:5]]
                            line_split[2:5] = [f"{coord:10.4f}" for coord in line_split[2:5]]
                        else:
                            line_split[column] = f"{line_split[column]:>10}"
                    lines[j] = ''.join(line_split) + '\n'
                    atom_idx += 1
                break
            if lines[i].startswith('@<TRIPOS>MOLECULE') and self._name_changed_:
                lines[i+1] = f"{self.mol.GetProp('_Name'):>7}\n"
        with open(out_file, 'w') as f:
            f.write("".join(lines))
        # import openbabel
        # obConversion = openbabel.OBConversion()
        # obConversion.SetInAndOutFormats("mol", "mol2")
        # obabel_mol = openbabel.OBMol()
        # obConversion.ReadString(obabel_mol, Chem.MolToMolBlock(self.mol))
        # obConversion.WriteFile(obabel_mol, out_file_name)
        # for idx, atom in enumerate(self.mol.GetAtoms()):
        #     obabel_mol.GetAtom(idx +1).SetPartialCharge(atom.GetDoubleProp('charge'))
        # obConversion.WriteFile(obabel_mol, 'taka')

def get_atom_index(file_path, H_atoms = True):
    """

    Parameters
    ----------
    file_path : TYPE, str
        DESCRIPTION. The file with the atoms, Could be an itp, gro or pdb.
    H_atoms : TYPE, optional
        DESCRIPTION. The default is True. :
            True:
                Return all the list of atoms.
            False:
                Just return the heavy atoms (non hydrogens).

    Returns
    -------
    atom_index : TYPE
        DESCRIPTION.

    """
    cwd = os.getcwd()
    tmpdir = tempfile.TemporaryDirectory()
    #Guessing the extension file
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"The file {file_path} doesn't exist or is not accessible")
    file_path = os.path.abspath(file_path)
    basename = os.path.basename(file_path)
    ext = basename.split('.')[-1].lower()
    os.chdir(tmpdir.name)
    atom_index = []
    if ext == 'itp':
        atoms_section = get_top_sections(file_path, dictionary=True)['atoms']
        for line in atoms_section:
            if not line.startswith(";") and not line == "\n":
                split = line.split()
                if H_atoms:
                    atom_index.append(int(split[0]))
                else:
                    if split[4].strip()[0] not in 'Hh': #This look for the first character of the atom name, if it is H or h, then is a hydrogen atom
                        atom_index.append(int(split[0]))


    elif ext in ['gro', 'pdb']:
        with open("ndx.opt", "w") as opt:
            opt.write("""
                      "System" group System;
                      "System-H" group System and ! name "H*" and ! name "h*"
                      """)
        run(f"""
                  echo 'q'| gmx make_ndx -f {file_path} -o tmp.ndx
                  gmx select -s {file_path} -sf ndx.opt -n tmp.ndx -on tmp
                  """)
        #deleting the line _f0_t0.000 in the file
        with open("tmp.ndx", "r") as index:
            data = index.read()
            data = data.replace("_f0_t0.000","")
        with open("tmp.ndx", "w") as index:
            index.write(data)

            ndx = NDX("tmp.ndx").data
            if H_atoms:
                atom_index = ndx["System"]
            else:
                atom_index = ndx["System-H"]


    # Returning to the original directory
    os.chdir(cwd)
    return atom_index
 # Both functions could be put together in a class NDX, that has one method, write, because always read. But this create the disatventage to incompatibiliti with ndx rectifivcation
 # En otras palabras que se pierde la ventaja de copiar un dictionario nen formatop NDX, se tendria que crear un NDX vacio (un archivo vacio , porque siempre lo va a leer) y asignar la data al objecto,

class NDX:
    def __init__(self, file):
        self.file = file
        self.data = {}
        self.read()


    def read(self):
        """
        Da bateo cuando el nombre del indice no tiene ningu atomo, dime algo!!!!
        Tengo que revisar en todos lo lugares que llamo a NDX, porque le voy a quitar que convierta a entero
        Y tengo que ver como hago mas general la idea de indices sin atomos


        Returns
        -------
        dictionary : TYPE: dict
            A dictionary with the group index in the form
            { group_name1: [atom_idndex1, atom_idndex2, ...]
            group_name2: [atom_idndex1, atom_idndex2, ...]
            ...
                }

        """
        with open(self.file,'r') as f:
            lines = f.readlines()

        positions = [i for (i, line) in enumerate(lines) if line.startswith('[')]

        for i in range(len(positions)):
            key = lines[positions[i]][2:-3]
            self.data[key] = []
            try:
                chunk_lines = lines[positions[i]+1:positions[i+1]]
            except:
                chunk_lines = lines[positions[i]+1:]
            for line in chunk_lines:
                self.data[key] += [int(value) for value in line.split()]


    def write(self, out_name = "index.ndx", backup = True):
        """
        Write a GROMACS index file

        Parameters
        ----------
        dictionary : dict
            The keys of the dictionary are the names of the groups.
            The values are the corresponded list of atoms
        out_name : TYPE, str
            DESCRIPTION. The default is "index.ndx".
            The name of the output file

        Returns
        -------
        None.

        """
        data = ""
        items_per_row = 15
        for key in self.data:
            data += f"[ {key} ]\n"
            i = 0

            while i < len(self.data[key]):
                string = [str(item) for item in self.data[key][i:i+items_per_row]]
                i += items_per_row
                data += " ".join(string) + "\n"

        if backup: backoff(out_name)
        out_file = open(out_name,'w')
        out_file.write(data)
        out_file.close()


def get_top_sections(topology, dictionary = False):
    """
    The flag dictionary is available because if in the topology there are two
    sections with the same name, then the last one is the only one that will
    be save.
    IN the case that

    Parameters
    ----------
    topology : path to a gromacs topology file
        DESCRIPTION.
    dictionary : TYPE, optional
        DESCRIPTION. The default is False.
        If True instead of a list return a dict with keywords the sections
        If not, a list with elements [key, info], ...

    Returns
    -------
    sections : list or dict depend on dictionary
        DESCRIPTION.
        The info of the topology file

    """
    with open(topology,"r") as top:
        top_lines = top.readlines()

    sections = [] #is not possible to use dict because there is some repeated sections, like dihedrals, I don't know why
    for i in range(len(top_lines)):
        if top_lines[i].startswith("["):

            key = ""
            for char in top_lines[i]:
                if char == "]":
                    break
                if char not in "[ ":
                    key += char
            key.strip()


            info = []
            for j in range(i + 1, len(top_lines)):
                if top_lines[j].startswith("["):
                    i = j - 1
                    break
                info.append(top_lines[j])

            sections.append([key.strip(), info])

    if dictionary: return dict(sections)
    return sections

def KbT(absolute_temperature):
    """Return the value of Kb*T in kJ/mol

    Args:
        absolute_temperature (float): The absolute temperature in kelvin
    """

    return absolute_temperature*CTE.Kb

def beta(absolute_temperature):
    return 1 / KbT(absolute_temperature)

def BoltzmannFactor(E,T):
    return np.exp(-E/KbT(T))

def get_force_cte(T, sigma):
    return KbT(T)/sigma**2 # Kj/(mol [sigma^2])


def get_sigma(T, force_constant):
    return np.sqrt(KbT(T)/force_constant)

#======================================================================================================
# Here I will introduce some nice tools for working with advances python modules for molecular dynamic
#======================================================================================================
def my_guess_types(names):
    """Guess the type of the atoms. It will simply remove the number of each string. Useful for the guessing of the ligands atoms

    Args:
        names (list): a list of atoms names
    """
    def f(array):
        for item in array:
            return ''.join(filter(str.isalpha, item))
    return np.vectorize(f)(names)

def array_to_bv(s):
    bv = ExplicitBitVect(len(s))
    on_bits = np.where(s == True)[0].tolist()
    bv.SetBitsFromList(on_bits)
    return bv

class classifier_ifp:
    def __init__(self, tpr, xtc, ligand_selection = 'resname LIA', protein_selection = 'protein', elapsed_steps = 0, load_fingerprint = False, ifp_file = 'ifp.pkl', ifp_atoms_file = 'ifp_atoms.pkl', n_clusters = 2) -> None:
        self.tpr = tpr
        self.xtc = xtc
        self.ligand_selection = ligand_selection
        self.protein_selection = protein_selection
        self.elapsed_steps = elapsed_steps
        self.load_fingerprint = load_fingerprint
        self.ifp_file = ifp_file
        self.ifp_atoms_file = ifp_atoms_file
        self.n_clusters = n_clusters
        self.parse()

    def __split(self, string):
        def f(s):
            if 'chain' in s:
                return s.split('_')[-1]
            elif "LI" in s:
                return s.split('_')[-1][-1]
            else:
                return ""
        return np.vectorize(f)(string)

    def cluster(self, n_clusters=2, max_iter=500, init='Cao', n_init=40, verbose=0, random_state=None, n_jobs=1):
        self.n_clusters = n_clusters
        ## Becasue all the data is True or False, every variable contribute with the same weight, so there is no need to rescale it.
        ## Initialize the model
        self.model = KModes(
            n_clusters=self.n_clusters,
            max_iter=max_iter,
            init=init,
            n_init=n_init,
            verbose=verbose,
            random_state=random_state,
            n_jobs=n_jobs
        )
        ## Run the model
        self.model.fit(self.ifp_df)
        # Group frames on clusters
        self.membership = {}
        for i in range(self.n_clusters):
            self.membership[i] = np.where(self.model.labels_ == i)[0]

    def parse(self):
        if self.load_fingerprint:
            self.ifp_df = pd.read_pickle(self.ifp_file).astype(np.uint8)
            try:
                self.ifp_df_atoms = pd.read_pickle(self.ifp_atoms_file)
            except:
                print('This object does not contains the atribute self.ifp_df_atoms')
        else:
            tmp_pdb = tempfile.NamedTemporaryFile(suffix='.pdb')
            run(f"export GMX_MAXBACKUP=-1;  gmx editconf -f {self.tpr} -o {tmp_pdb.name}")


            ref = mda.Universe(tmp_pdb.name) # This has the correct numeration for the residues
            u = mda.Universe(self.tpr, self.xtc)#, guess_bonds=True

            guessed_elements = guess_types(my_guess_types(u.atoms.names))
            u.add_TopologyAttr('names', guessed_elements)
            u.add_TopologyAttr('elements', guessed_elements)
            u.add_TopologyAttr('types', guessed_elements)
            u.add_TopologyAttr('resid', ref.residues.resids)
            u.add_TopologyAttr('chainID', self.__split(u.atoms.segids)) # Puede pensar en algo mas general con la funcion que tengo dentro del pdb para adivinar las cadenas

            self.ligand = u.atoms.select_atoms(self.ligand_selection)
            self.protein = u.atoms.select_atoms(self.protein_selection)

            # Initialize fingerprint
            ifp = plf.Fingerprint()
            # run on a slice of frames from begining to end with a step of 10
            if self.elapsed_steps:
                ifp.run(u.trajectory[::self.elapsed_steps], self.ligand, self.protein)
            else:
                ifp.run(u.trajectory, self.ligand, self.protein)

            self.ifp_df = ifp.to_dataframe()
            self.ifp_df_atoms = ifp.to_dataframe(return_atoms = True)

            # Saving the fingerprint
            self.ifp_df.to_pickle(self.ifp_file)
            # Also the ifp_atoms will be exported
            self.ifp_df_atoms.to_pickle(self.ifp_atoms_file)
            tmp_pdb.close()

        # Run cluster with the default values nad number of cluster
        self.cluster(n_clusters=self.n_clusters)

    def get_closers(self, user_membership = None):
        """Return the closest frames for each cluster in self.membership or user_memebership if provided.

        Args:
            user_membership (dict, optional): The memebership, the keys will be the cluster identifier and the values the frames taht belongs to each cluster. Defaults to None.

        Returns:
            _type_: _description_
        """
        bvs_centroids = [array_to_bv(centroids) for centroids in self.model.cluster_centroids_]
        bvs_frames = plf.to_bitvectors(self.ifp_df)

        if user_membership:
            membership = user_membership
        else:
            membership = self.membership

        closest_frames = {}
        for cluster in membership:
            bvs_selected_frames = [bvs_frames[frame] for frame in membership[cluster]]
            # If a cluster only have one member
            if len(bvs_selected_frames) == 1:
                # BulkTanimotoSimilarity needs more than one fingerprint, this is just a trick, and probably alwasy be avoid.
                bvs_selected_frames += bvs_selected_frames
            similarities = DataStructs.BulkTanimotoSimilarity(bvs_centroids[cluster], bvs_selected_frames)
            index = similarities.index(max(similarities))
            closest_frames[cluster] = membership[cluster][index]
        return closest_frames

class classifier_dist_contacts:
    def __init__(self, tpr, xtc, chains = ('A', 'B', 'C', 'D', 'E'), elapsed_steps = 0, n_clusters = 2) -> None:
        self.tpr = tpr
        self.xtc = xtc
        self.chains = chains
        self.elapsed_steps = elapsed_steps
        self.n_clusters = n_clusters
        self.data = pd.DataFrame()
        self.parse()

    def cluster(self, n_clusters=2, init='k-means++', n_init=50, max_iter=500, tol=0.0001, verbose=0, random_state=None, copy_x=True, algorithm='auto'):
        features = ['AveNumResi', 'AveDistZ']
        self.n_clusters = n_clusters
        pipe = Pipeline(
            [
                ("scaler", MinMaxScaler()),
                ('kmeans',
                    KMeans(
                        n_clusters=self.n_clusters,
                        init=init,
                        n_init=n_init,
                        max_iter=max_iter,
                        tol=tol,
                        verbose=verbose,
                        random_state=random_state,
                        copy_x=copy_x,
                        algorithm=algorithm
                    )
                )
            ]
        )
        pipe.fit(self.data[features])
        self.model = pipe['kmeans']
        #Think about if implement all the method in a athomatic way to improve the selction of the clusters
        # Group frames on clusters
        self.membership = OrderedDict()
        for i in range(self.n_clusters):
            self.membership[i] = np.where(self.model.labels_ == i)[0]

    def parse(self):
        u = mda.Universe(self.tpr, self.xtc)

        close_contacts_frame_0 = {chain:u.atoms.select_atoms(f'(backbone and segid *Protein_chain_{chain}) and around 4 resname LI{chain}') for chain in self.chains}

        if self.elapsed_steps:
            iterator = u.trajectory[::self.elapsed_steps]
            Frame = range(0,u.trajectory.n_frames,self.elapsed_steps)
        else:
            iterator = u.trajectory
            Frame = range(u.trajectory.n_frames)
        AveNumResi = []
        AveDistZ = []
        for ts in tqdm.tqdm(iterator):
            NumResi = []
            DistZ = []
            for chain in self.chains:
                NumResi.append(len(u.atoms.select_atoms(f'(protein and segid *Protein_chain_{chain}) and around 3 resname LI{chain}').residues))
                DistZ.append(np.abs(close_contacts_frame_0[chain].center_of_mass(pbc=True, compound='group')[-1] - u.atoms.select_atoms(f'resname LI{chain}').center_of_mass(pbc=True, compound='group')[-1]))
            AveNumResi.append(np.mean(NumResi))
            AveDistZ.append(np.mean(DistZ))
        self.data = pd.DataFrame(
            {
                'Frame': Frame,
                'AveNumResi': AveNumResi,
                'AveDistZ': AveDistZ,
            }
        )
        # Run cluster with the default values nad number of cluster
        self.cluster(n_clusters=self.n_clusters)

    def get_closers(self):
         # https://stackoverflow.com/questions/26795535/output-50-samples-closest-to-each-cluster-center-using-scikit-learn-k-means-libr
        # This is like on the paper was done, but I need to do the scalling first
        scale = MinMaxScaler()
        scale.fit(self.data[['AveNumResi', 'AveDistZ']])
        closest_frames, distances = vq(self.model.cluster_centers_, scale.transform(self.data[['AveNumResi', 'AveDistZ']]))
        return closest_frames

if __name__ == '__main__':

    pass
    print(CTE.Kb)
    #k = aovec([(0,0,1), (0,1,0), (1,0,0)])
    #deg = 2.5
    #rad = np.pi*deg/180
    #n = 0.05
    #print(get_force_cte(298.15, n))
    #print(get_sigma(303.15, 500))
    # now = datetime.datetime.now()
    #ndx = NDX("/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_J/7e27/sys_MMV007839_Cell_891_SP_param/index.ndx")
    #print(ndx.data)
    # check = CHECK("md.log", "md.lis")
    # print(check.performance)
    # print(check.estimated_end_datetime())
    # print(check.daysdiff(now))
    #print(os.listdir())
    #cp("test.py", "h/k.py")
    #print(list_if_file("./"))q
    #print(os.path.isfile("./*py"))
    # tpr = '/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_M_annealing/7e27/sys_MMV007839_Cell_891_SP_param/pull.tpr'
    # xtc = '/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_M_annealing/7e27/sys_MMV007839_Cell_891_SP_param/pull.xtc'
    # fingerprint_files = [
    #     '/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_M_annealing/7e27/sys_MMV007839_Cell_891_SP_param/ifp.pkl',
    #     '/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_M_annealing/7e27/sys_MMV007839_Cell_891_SP_param/ifp_atoms.pkl'
    # ]
    # ligand = 'LIA'

    # a = classifier_dist_contacts(tpr, xtc, elapsed_steps=200)
    # #print(a.data)
    # #print(a)
    # #a = classifier_ifp(tpr, xtc, ligand_selection='resname LIA', elapsed_steps = 0, load_fingerprint = True, ifp_file = fingerprint_files[0], n_clusters = 10)
    # print(a.model.cluster_centers_.shape)
    # a.cluster(n_clusters=3)

    # print(a.model.cluster_centers_)
    # print(a.get_closers())
    # #print(a.get_closers())
    # #print(a.labels_)
    # print(KbT(300))
    # print(get_force_cte(298,1))
    # d0 = xvg.XVG('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_N_ST/7e27/sys_MMV007839_Cell_891_SP_param/windows/coord0_selected.xvg').data
    # d1 = xvg.XVG('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_N_ST/7e27/sys_MMV007839_Cell_891_SP_param/windows/coord1_selected.xvg').data
    # d2 = xvg.XVG('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_N_ST/7e27/sys_MMV007839_Cell_891_SP_param/windows/coord2_selected.xvg').data
    # d3 = xvg.XVG('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_N_ST/7e27/sys_MMV007839_Cell_891_SP_param/windows/coord3_selected.xvg').data
    # d4 = xvg.XVG('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_N_ST/7e27/sys_MMV007839_Cell_891_SP_param/windows/coord4_selected.xvg').data
    # d5 = xvg.XVG('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_N_ST/7e27/sys_MMV007839_Cell_891_SP_param/windows/coord5_selected.xvg').data
    # data = [d0, d1, d2, d3, d4, d5]
    # AlignPoint = 4.5
    # NumbPoints = 300
    # print(RSS(data,AlignPoint, NumbPoints = NumbPoints))