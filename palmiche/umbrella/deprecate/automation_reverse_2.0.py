#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Fri Jun 11 19:31:35 2021
Author        : Alejandro Martínez León
Mail          : [alejandro.martinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""
"""
The correct way to use this script is from the above directory and 
htmd, to activate the environment
nohup ./auto_umbrella.py >output.txt 2>&1 &
"""




from palmiche.umbrella import templates, COM_distance, conf_windows
from palmiche.utils import tools, mdp
import os, glob, time, datetime, tempfile, argparse, re
import numpy as np


def get_conf_from_type(conf_type, path, receptor, ligand):
    """This is used for dict_gen. It is just to select whta will be the initial configuration for the reverse pulling: the last frame of the pull (pull) or the last equilibrated window (window).

    Args:
        conf_type ([type]): [description]
        path ([type]): [description]
        receptor ([type]): [description]
        ligand ([type]): [description]

    Returns:
        [type]: [description]
    """
    if conf_type == 'pull':
        conf = os.path.join(path, receptor, ligand, 'pull.gro')
    elif conf_type == 'window':
        pattern = re.compile("^\d{5}$")
        windows_path = os.path.join(path, receptor, ligand, 'windows')
        maximum = max([os.path.join(windows_path, window) for window in tools.list_if_dir(windows_path) if pattern.match(window)])
        conf = os.path.join(maximum, 'npt.gro')
    return conf

def dict_gen(path = "./umbrella", FF_protein = "amber99sb-star-ildn.ff",
             FF_lipid="Slipids_2020.ff", FF_ligand="GAFF2.ff", topol="topol_heavyH.top",
             toppar="toppar", conf_type = "window", ndx = 'index.ndx'):
    """
    Mejorar descripcion
    THos function demand to have a tree shape:
        path/receptor1/ligand1 and so on
        path/receptor2/ligand456 and so on and with the files or directories
        defined in the function inside the last directory (ligand

    Parameters
    ----------
    path : TYPE, optional
        DESCRIPTION. The default is "./HMR_2_5".
    FF_protein : TYPE, optional
        DESCRIPTION. The default is "amber99sb-star-ildn.ff".
    FF_lipid : TYPE, optional
        DESCRIPTION. The default is "Slipids_2020.ff".
    FF_ligand : TYPE, optional
        DESCRIPTION. The default is "GAFF2.ff".
    topol : TYPE, optional
        DESCRIPTION. The default is "topol_heavyH.top".
    toppar : TYPE, optional
        DESCRIPTION. The default is "toppar".
    conf : TYPE, optional
        DESCRIPTION. The default is "confout.gro".

    Returns
    -------
    results : TYPE
        DESCRIPTION.

    """


    
    path = os.path.abspath(path)

    dict_paths = {}
    for receptor in tools.list_if_dir(path):
        # In this way if the receptor where teasted to different ligands  (een so, they must be in the deffinition of the ligands_dict) this will tak it intop account
        for ligand in tools.list_if_dir(os.path.join(path, receptor)):
            dict_paths.setdefault(receptor, {})[ligand] = {"FF_protein": os.path.join(path, receptor, ligand, FF_protein),
                                                           "FF_lipid": os.path.join(path, receptor, ligand, FF_lipid),
                                                           "FF_ligand": os.path.join(path, receptor, ligand, FF_ligand),
                                                           "topol": os.path.join(path, receptor, ligand, topol),
                                                           "toppar": os.path.join(path, receptor, ligand, toppar),
                                                           "conf": get_conf_from_type(conf_type, path, receptor, ligand),
                                                           "ndx": os.path.join(path, receptor, ligand, ndx),
                                                           }
    return dict_paths


def main(input_path_dict,
         output_path,
         dt=0.004,
         simulation_time={'pulling': 100000,
                          'equilibration': 2000, 'production': 20000},#ps
         number_frame={'pulling': 2000,
                       'equilibration': 2, 'production': 2000},#ps
         temperature = 303.15,
         annealing = True,
         pull_coord_vec = [0.0, 0.0, -1.0],
         flat_bottom_init=0.5,  # This is the radium
         flat_bottom_k=400,
         orient_restrain_init=20,  # 0.349, #; 20 degrees 20/180*3.141
         orient_restrain_k=2.5,  # #; in kJ/mol/deg^2, leads to sigma = 1 degrees k = KbT/sigma^2
         # This means that the vector will be calculate in suing the initial orientation fo the ligands. if a vector is provied, for example (0,0,-1), then this will be used
         orient_restrain_vec="calculate",
         pull_distance=4.5,
         # Here I will identify what are the different intervals on the windows.
         pull_middle_points=[],
         # Here I will deffine how many pulling force constant I will use, if there is only on number then the same applied for every thong, I could even add the pulling force constant for the pulling itself and the one for the windows
         pull_force_constants=[1500],
         # THis variable also will have as many as different intervals are specified, in case only one, tha same applied for all the windows
         window_widths=[0.1],
         chains=['A', 'B', 'C', 'D', 'E'],
         GROMACS_version="2021.1",
         cpu=12,
         COM_dist_cpu=12,
         machine="smaug",
         partition="deflt",
         exclude_nodes="",
         pull_code=True,
         window_code=True):
    """
    The three values of simulation_time, number_frame, flatt_bottom_r and flatt_bottom_k
    are because we need to generate
    three different mdp:
    1- Pulling
    2-NPT equilibration of each window
    3-MD of each window

    Parameters
    ----------
    input_path_dict : TYPE, dictionary of path
        DESCRIPTION. This dictionary is refereed to the one that you get from:
            dict_gen()['receptor_name']['ligand_name']
    output_path : TYPE, output path 
        DESCRIPTION. Destination      
    dt : TYPE, optional
        DESCRIPTION. The default is 0.004. [ps]
    simulation_time : TYPE, optional, dict
        DESCRIPTION. The default is {'pulling':100,'equilibration':1,'production':10}. [ns]
    number_frame : TYPE, optional
        DESCRIPTION. The default is {'pulling':2000,'equilibration':2,'production':2000}.
    flat_bottom_r : TYPE, optional
        DESCRIPTION. The default is {'pulling':0.3,'equilibration':0.4,'production':0.5. [nm]
    flat_bottom_k  : TYPE, optional
        DESCRIPTION. The default is {'pulling':600,'equilibration':500,'production':400}. [kJ/mol]
    pull_distance : TYPE, optional
        DESCRIPTION. The default is 3.5. [nm]
    pull_middle_points : TYPE, list optional
        DESCRIPTION. The default is []. Here we deffine the middle points for the windows, In case that any, you must provied an empty list
    pull_force_constants  : TYPE, list optional
        DESCRIPTION. The default is [100000] kJ/mol. Here we deffine the force constants. You could provied:
            1-) one: will be used for all simulations).
            2-) two: the first one for the pulling and the second one for all the windows
            3-) n: In this case you need to specified the force constant for the pulling (the first o the list), and all for each window.
    window_widths : TYPE, list optional
        DESCRIPTION. The default is [0.1]. The width for the windows. This could have:
            1-) Only one value, and will be applied for all the windows
            2-) As many as interval were deffined with pull_middle_points
    chains : TYPE, optional
        DESCRIPTION. The default is ["A", "B","C", "D", "E"]
    itp_ligands : TYPE, optional, str (glob expression) or list
        DESCRIPTION. The default is "LI['A', 'B', 'C', 'D', 'E'].itp". 
        This is the name of the ligands itp files in the toppar directory.
        In this case the string will be interpreted as a possible glob expresion
        but if a list is provided, then will look for the items explicitely decleared.
        E.g ['LIA', 'LIB','LIC', etc...]
    GROMACS_version : TYPE, optional
        DESCRIPTION. The default is "2021.1".
    cpu : TYPE, optional
        DESCRIPTION. The default is 12.
        Number of logical cpu to use for gmx mdrun
    COM_dist_cpu : TYPE, optional
        DESCRIPTION. The default is 12.
        Number of logical cpu to use for COM_distance.main()
    slurm : TYPE, optional
        DESCRIPTION. The default is smaug.
        smaug, gwdg, frontend
    Tengo que completar la documentacion
    Returns
    -------
    None.

    """
    # Here I am see how many different intervals are in the windows
    # Selecting the correct set of parameters for pull_force_constants and window_widths
    pull_middle_points = sorted(pull_middle_points)
    print(
        f"The following window intervals will be constructed:\nmin---{'---'.join([str(point) for point in pull_middle_points])}---max\n")
    if len(pull_force_constants) == 1:
        pull_force_constants = (
            len(pull_middle_points) + 2) * pull_force_constants
        print(
            f"The force constant {pull_force_constants[0]} will be used for all the simulations.\n")
    elif len(pull_force_constants) == 2:
        pull_force_constants = [pull_force_constants[0]] + [pull_force_constants[1]
                                                            for i in range(len(pull_middle_points) + 1)]
        print(
            f"The force constant {pull_force_constants[0]} will be used for the pulling and {pull_force_constants[1]} for the production simulation  of all windows.\n")
    elif len(pull_force_constants) != len(pull_middle_points) + 2:
        raise SyntaxError(f"""
    pull_force_constants could have:
            1-) Only one value. In this case the same force constant will be applied to all the simulations.
            2-) Two values. In this case the first one will be used for the pulling and the second one for all the windows.
            3-) len(pull_middle_points) + 2 . In this case the force constant is specified for each simulation. The first one for the pulling, the second one for the first interval of windows, the third one for the second one and so on. In the case that any middle points were specified pull_force_constants could only have one or two values.
    However pull_middle_points = '{pull_middle_points}' and  pull_force_constants = '{pull_force_constants}'.
    """)
    else:
        print(
            f"The force constant {pull_force_constants} will be used. The first one is for the pulling simulation, the rest one for the windows.\n")

    if len(window_widths) == 1:
        window_widths = (len(pull_middle_points) + 1) * window_widths
        print(
            f"The window widths {window_widths[0]} will be used for all windows")
    elif len(pull_middle_points) + 1 != len(window_widths):
        raise SyntaxError(f"""
    window_widths could have:
            1-) Only one value. In this case the same window interval will be applied.
            2-) len(pull_middle_points) + 1 . In this case the window interval is specified for each simulation. The first one for first interval of windows, the second one for the second one window and so on. In the case that any middle points were specified window_widths could only have one value.
    However pull_middle_points = '{pull_middle_points}' and  window_widths = '{window_widths}'.
    """)
    else:
        print(f"The window widths {window_widths} will be used.\n")

    # In case of more partition in the cluster you need to provided to the dictionary

    maximum_days = {"deflt": 2,
                    "long": 5}

    cwd = os.getcwd()
    start_datetime = datetime.datetime.now()
    output_path = os.path.abspath(output_path)

    # Generating the ligands name whit the chain information
    ligands = [f"LI{chain}" for chain in chains]
    # Generating the tc_grps names
    tc_grps = ['SOLU', 'MEMB', 'SOLV'] + ligands
    # The general options for the mdp
    mdp_options = {
        "dt": dt,
        "tc_grps": ' '.join(tc_grps),
        "tau_t": ' '.join(len(tc_grps)*['1.0']),
        "ref_t": ' '.join(len(tc_grps)*[str(temperature)])
    }



    window_path = os.path.join(output_path, 'windows')
    new_path_dict = {}

    for key in input_path_dict:
        if key != 'conf':
            new_path_dict[key] = os.path.join(output_path, os.path.basename(input_path_dict[key]))


    [tools.makedirs(item) for item in [output_path, window_path]]

    # The two if for pull_code and for window_code is in order to flexibly the running of the script
    if pull_code:
        # Ordering the files, now new_path_dict is useful.
        [tools.cp(input_path_dict[key], output_path, r=True)  for key in input_path_dict if key != 'conf']

        # Changing the directory to
        os.chdir(output_path)

        # Creating the average oriented vector. This vector is used for the orientation restrain of the ligands. This is done is order to have
        # more consisten the simulation and the stadistic better for the 5 PMF, because if not he ligands have different orientation, and that means that the PMFs are not comparables.
        # I will calculate only fi the vector provi
        try:  # I need to use try, because if a list was provided then, will raise an error, because the list doesn't have strip() method
            if orient_restrain_vec.strip().lower() == "calculate":
                # first I need to create a temporal tpr file
                mdptmp = tempfile.NamedTemporaryFile(suffix=".mdp")
                tprtmp = tempfile.NamedTemporaryFile(suffix=".tpr")
                tools.run(
                    f"export GMX_MAXBACKUP=-1; gmx grompp -f {mdptmp.name} -c {os.path.basename(input_path_dict['conf'])} -r {os.path.basename(input_path_dict['conf'])} -p {os.path.basename(input_path_dict['topol'])} -n {input_path_dict['ndx']} -o {tprtmp.name}")

                # Here I am calculating all the vectors for each ligand
                vectors = [tools.get_vec_COM(os.path.basename(input_path_dict['conf']), tprtmp.name,
                                             input_path_dict['ndx'], f"{ligand}_GROUP1", f"{ligand}_GROUP2") for ligand in ligands]
                # Get the average oriented vector, rounded with 3 decimals
                orient_restrain_vec = tools.aovec(vectors, round=3)
        except:
            pass
        # Generating the mdp using the templates module, here we are constructing the pulling.
        MDP_pulling = mdp.MDP(
            **mdp_options,
            define = f"-DPOSRES -DPOSRES_FC_BB=1000.0 -DPOSRES_FC_SC=100.0 -DPOSRES_FC_LIPID=1000.0 -DDIHRES -DDIHRES_FC=1000.0 -DPOSRES_LIG=0.0",
            time = simulation_time['pulling'],
            xtc_numb_frame = number_frame['pulling'])
        MDP_pulling.pull(
            pull_distance,
            ligands,
            # Options for the cylinder and orientation restraints
            flat_bottom_init=flat_bottom_init,
            flat_bottom_k=flat_bottom_k,
            orient_restrain_init=orient_restrain_init,
            orient_restrain_k=orient_restrain_k,
            orient_restrain_vec=orient_restrain_vec,
            # General options for the actual pulling
            pull_coord1_k=pull_force_constants[0],
            pull_coord1_vec=" ".join([str(xi) for xi in pull_coord_vec]),
            pull_coord1_init=pull_distance,
            pull_coord1_rate= - pull_distance / (simulation_time['pulling']), # I need to give the negative rate, this is calculated automatically by the module by I need to provied the negative sign, for that reason I calculated here
            )
        MDP_pulling.write(os.path.join(output_path, 'pull.mdp'))

        # Creating the tpr file
        tools.run(f"gmx grompp -f pull.mdp -c {input_path_dict['conf']} -r {input_path_dict['conf']} -p {input_path_dict['topol']} -n {input_path_dict['ndx']} -o pull.tpr")

        with open(os.path.join(output_path, "job.sh"), "w") as file:
            file.write(templates.GROMACS_PULL_JOB(
                job_name=f"pull",
                tpr_name="pull",
                partition=partition,
                cpus_per_task=cpu,
                nodes=1,
                nice=0,
                gpus=1,
                version=GROMACS_version,
                exclude=exclude_nodes,
                gromacs_cmd_add="-update gpu",
                machine=machine))
        # Executing the job
        if machine != "frontend":
            JOBID_tmp_file = tempfile.NamedTemporaryFile()
            tools.run(f"sbatch --parsable job.sh>{JOBID_tmp_file.name}")
            start_tmp_datetime = datetime.datetime.now()
            with open(JOBID_tmp_file.name, "r") as f:
                JOBID = int(f.readline())
            while JOBID in tools.checkrun():
                time.sleep(20)
            # In case tha the job ended. Now we will check if the job ended:
            #   1- If check.performance is True means that the job ended (None is the value in case that the simulation did not end)
            #   2- If not but the time estimated is longer than the one offered by the selected partition, run as many time needed
            #           (this could give natural errors in the case that the time is longer but the simulation finished for other reason)
            #   3- When all the possible simulation finish, if still not ended correctly raise an error, It will not go to the next steps
            # I use md.lis because is the one used by default in jobscript_slurm.
            check = tools.CHECK("pull.log", "md.lis")
            if not check.performance:
                daysdiff = check.daysdiff(start_tmp_datetime)
                if daysdiff > maximum_days[partition]:
                    # This will give me how many times I need to run again the simulation in order to complet it (-1 because I already run one time)
                    for i in range(int((daysdiff / maximum_days[partition]) - 1)):
                        # Run again the job
                        JOBID_tmp_file = tempfile.NamedTemporaryFile()
                        tools.run(
                            f"sbatch --parsable job.sh>{JOBID_tmp_file.name}")
                        with open(JOBID_tmp_file.name, "r") as f:
                            JOBID = int(f.readline())
                        while JOBID in tools.checkrun():
                            time.sleep(20)
                        check = tools.CHECK("pull.log", "md.lis")
                        if check.performance:
                            break
                else:
                    raise RuntimeError("The pull simulation ended with error.")
            # Check for the last time if the simulation ended correctly, if not, it is not due to time issues. capabilities
            # I use md.lis because is the one used by default in jobscript_slurm.
            check = tools.CHECK("pull.log", "md.lis")
            if not check.performance:
                raise RuntimeError("The pull simulation ended with error.")
        else:
            tools.run(f"bash job.sh")

    if window_code:

        # Changing dir to windows and measurement of distances
        os.chdir(window_path)
        grouplist = [[ligand, ligand+"_CLOSE_AA"] for ligand in ligands]
        # This will output a file called mean.dat, that has the format needed for pulling_processing.windows_creator. The distance is the average distance of the grouplist
        COM_distance.main(grouplist, cpu=COM_dist_cpu, ndx=new_path_dict['ndx'], tpr=os.path.join(output_path,'pull.tpr'),
                          xtc=os.path.join(output_path, 'pull.xtc'), prefix='conf', out='summary_distances.dat')

        # ==============================================================================
        # The processing and the launch of the windows

        # ==============================================================================

        # _______checking__________
        #check_cont = input("Now the configuration of the windows will take place. Will you continue?[y/n]:").lower()[0]
        #check_cont = "y"
        # _________________________
        # if check_cont == "y":
        # This part of the code is for the creation of the windows with their respective pulling forces and interval widths.
        windows = {}
        if pull_middle_points:
            maximum_windows_ID = 0
            for (i, width) in enumerate(window_widths):
                if i == 0:
                    tmp_windows = conf_windows.main(
                        "mean.dat", width, txt_summary=None,  maximum=pull_middle_points[i])
                    for key in tmp_windows:
                        windows[key] = tmp_windows[key]
                        windows[key]['interval'] = f"min-{pull_middle_points[i]}"
                        windows[key]['width'] = width
                        # Because the first force_constant is for the pulling
                        windows[key]['force_constant'] = pull_force_constants[i+1]
                elif i == len(window_widths) - 1:
                    # Here is added by default the sampled distance to the maximum value
                    tmp_windows = conf_windows.main(
                        "mean.dat", width, txt_summary=None,  minimum=pull_middle_points[i - 1])

                    for key in tmp_windows:
                        # Only if the frame is not already in the dict, I will add the window
                        if tmp_windows[key]['frame'] not in [windows[key2]['frame'] for key2 in windows]:
                            windows[key + maximum_windows_ID] = tmp_windows[key]
                            windows[key +
                                    maximum_windows_ID]['interval'] = f"{pull_middle_points[i - 1]}-max"
                            windows[key + maximum_windows_ID]['width'] = width
                            # Because the first force_constant is for the pulling
                            windows[key + maximum_windows_ID]['force_constant'] = pull_force_constants[i+1]
                        else:
                            # If the frame is repeated I need to substract the len of the tmp_window
                            maximum_windows_ID -= 1
                else:
                    tmp_windows = conf_windows.main(
                        "mean.dat", width, txt_summary=None,  minimum=pull_middle_points[i - 1], maximum=pull_middle_points[i])

                    for key in tmp_windows:
                        # Only if the frame is not already in the dict, I will add the window
                        if tmp_windows[key]['frame'] not in [windows[key2]['frame'] for key2 in windows]:
                            windows[key + maximum_windows_ID] = tmp_windows[key]
                            windows[key +
                                    maximum_windows_ID]['interval'] = f"{pull_middle_points[i - 1]}-{pull_middle_points[i]}"
                            windows[key + maximum_windows_ID]['width'] = width
                            # Because the first force_constant is for the pulling
                            windows[key + maximum_windows_ID]['force_constant'] = pull_force_constants[i+1]
                        else:
                            # If the frame is repeated I need to substract the len of the tmp_window
                            maximum_windows_ID -= 1
                maximum_windows_ID += len(tmp_windows)
        else:
            windows = conf_windows.main(
                "mean.dat", window_widths[0], txt_summary=None)
            for key in windows:
                windows[key]['interval'] = "min-max"
                windows[key]['width'] = window_widths[0]
                # Because the first force_constant is for the pulling
                windows[key]['force_constant'] = pull_force_constants[1]

        out_string = f"{'window_ID':<15}{'frame':<10}{'dist':<10}{'delta_dist':<15}{'interval':<15}{'width':<10}{'force_constant':<10}\n"
        for key in windows:
            out_string += f"{key:<15}{windows[key]['frame']:<10}{windows[key]['dist']:<10.3}{windows[key]['delta_dist']:<15}{windows[key]['interval']:<15}{windows[key]['width']:<10.3}{windows[key]['force_constant']:<10}\n"
        with open("conf_windows_summary.txt", 'w') as o:
            o.write(out_string)

        tools.makedirs("split_xtc")
        tpr_files = open(os.path.join(window_path, "tpr_files.dat"), "w")
        pullf_files = open(os.path.join(window_path, "pullf_files.dat"), "w")
        pullx_files = open(os.path.join(window_path, "pullx_files.dat"), "w")

        for conf in glob.glob("conf[0-9]*"):

            # Here we are taken the frame and see which is the window ID of this frame
            frame = int(conf.split(".")[0][4:])
            ID = None
            for key in windows:
                if windows[key]['frame'] == frame:
                    ID = key
                    break

            if ID:
                ID_dir_basename = str(ID).zfill(5)
                tools.makedirs(ID_dir_basename)
                ID_dir_abspath = os.path.abspath(ID_dir_basename)
                tools.mv(conf, ID_dir_abspath)

                [tools.cp(new_path_dict[key], ID_dir_abspath, r=True) for key in new_path_dict]

                # Generating the mdp files.

                # First we define the general parameters for the NPT equilibration
                # ******************************************************************
                MDP_equilibration = mdp.MDP(
                    **mdp_options,
                    define=f"-DPOSRES -DPOSRES_FC_BB=100.0 -DPOSRES_FC_SC=10.0 -DPOSRES_FC_LIPID=100.0 -DDIHRES -DDIHRES_FC=100.0 -DPOSRES_LIG=0.0",
                    time=simulation_time['equilibration'],
                    xtc_numb_frame=number_frame['equilibration'])
                MDP_equilibration.pull(
                    pull_distance,
                    ligands,
                    # Options for the cylinder and orientation restraints
                    flat_bottom_init=flat_bottom_init,
                    flat_bottom_k=flat_bottom_k,
                    orient_restrain_init=orient_restrain_init,
                    orient_restrain_k=orient_restrain_k,
                    orient_restrain_vec=orient_restrain_vec,
                    # General options for the actual pulling
                    pull_coord1_k=windows[ID]['force_constant'],
                    pull_coord1_vec=" ".join([str(xi) for xi in pull_coord_vec]),
                    pull_coord1_init=windows[ID]['dist'],
                    pull_coord1_rate= 0.0,)

                # Adding the section for the annealing if required. 
                if annealing:
                    # This is a fix method
                    # Increase in 50 K the temperature for the ligands and the protein
                    # during the initial 25 % of the time, then the continued with a cooling (50%)
                    # And finally set constant the temperature to the initial one
                    tc_grps = np.array(tc_grps)
                    GroupBool = (tc_grps != 'MEMB') & (tc_grps != 'SOLV')                 
                    temp_list = np.arange(temperature, temperature + 52, 2)
                    MDP_equilibration.annealing2(GroupBool, temp_list, heating_frac=0.25, cooling_frac=0.50)

                MDP_equilibration.write(os.path.join(ID_dir_abspath, 'npt.mdp'))

                # Second we define the parameters for the umbrella production
                # ******************************************************************
                MDP_production = mdp.MDP(
                    **mdp_options,
                    define=f"-DPOSRES -DPOSRES_FC_BB=0.0 -DPOSRES_FC_SC=0.0 -DPOSRES_FC_LIPID=0.0 -DDIHRES -DDIHRES_FC=0.0 -DPOSRES_LIG=0.0",
                    time=simulation_time['production'],
                    xtc_numb_frame=number_frame['production'])
                MDP_production.pull(
                    pull_distance,
                    ligands,
                    # Options for the cylinder and orientation restraints
                    flat_bottom_init=flat_bottom_init,
                    flat_bottom_k=flat_bottom_k,
                    orient_restrain_init=orient_restrain_init,
                    orient_restrain_k=orient_restrain_k,
                    orient_restrain_vec=orient_restrain_vec,
                    # General options for the actual pulling
                    pull_coord1_k=windows[ID]['force_constant'],
                    pull_coord1_vec=" ".join([str(xi) for xi in pull_coord_vec]),
                    pull_coord1_init=windows[ID]['dist'],
                    pull_coord1_rate= 0.0,)
                
                MDP_production.write(os.path.join(ID_dir_abspath, 'production.mdp'))

                with open(os.path.join(ID_dir_abspath, "job.sh"), "w") as file:
                    file.write(templates.GROMACS_WINDOW_JOB(
                        job_name=f"window_{ID}",
                        npt_mdp_name="npt",
                        umbrella_mdp_name="production",
                        conf_file=conf,
                        posre_file=conf,
                        topol=os.path.basename(new_path_dict['topol']),
                        index=os.path.basename(new_path_dict['ndx']),
                        partition=partition,
                        cpus_per_task=cpu,
                        nodes=1,
                        nice=0,
                        gpus=1,
                        version=GROMACS_version,
                        exclude=exclude_nodes,
                        gromacs_cmd_add="-update gpu",
                        machine=machine))

                tpr_files.write(f"./{ID_dir_basename}/production.tpr\n")
                pullf_files.write(
                    f"./{ID_dir_basename}/production_pullf.xvg\n")
                pullx_files.write(
                    f"./{ID_dir_basename}/production_pullx.xvg\n")

            else:
                tools.mv(conf, "split_xtc")

        tpr_files.close()
        pullf_files.close()
        pullx_files.close()

        number_windows = len(windows)
        # THis case I use 3 because I use 2 pull restrain peer coord. orientation and loss of channel entrance
        number_coords = 3*len(ligands)
        select_file_names = []
        # THis assumes that the first coordinates are the one for the pulling and not the restrains coordinates.
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

        if machine != "frontend":
            JOBIDs = tools.job_launch()
            while len(list(set(JOBIDs).difference(tools.checkrun()))) != len(JOBIDs):
                time.sleep(20)
        else:
            JOBIDs = tools.job_launch(shell='bash')  # This is just for testing



        #wham_general = ["gmx wham -it tpr_files.dat -if pullf_files.dat -o all_coords_selected -hist hist_all_coords_selected -unit kJ"]
        wham_individual_coords = [f"gmx wham -ac -temp {temperature} -zprof0 {pull_distance} -bins 300 -unit kJ -nBootstrap 200 -bs-method hist -is {s} -ix pullx_files.dat -it tpr_files.dat -o {s.split('.')[0]} -hist hist_{s.split('.')[0]} -oiact iact_{s.split('.')[0]}.xvg -bsres bsResult_{s.split('.')[0]}.xvg -bsprof bsProfs_{s.split('.')[0]}.xvg" for s in select_file_names]
        for cmd in wham_individual_coords:  # wham_general + wham_individual_coords:
            tools.run(cmd)

    os.chdir(cwd)
    end_datetime = datetime.datetime.now()
    print(
        f"The function automation.main() finished. The script start at {start_datetime} and ended at {end_datetime}; the execution time was: {round((end_datetime - start_datetime).total_seconds() / 86400, 3)} days!")


if __name__ == '__main__':
    """This is just the main function of the script executing top_parser, atom_block_modifier
    and writer functions, see __doc__ for details"""

    #@@@@@@@@@@@@Esta parte tengo que mejorarla con los argumentos que se le pasan.

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--input_path_dict',
                        help="String-like-dict. This dictionary is refereed to the one that you get from: dict_gen()['receptor_name']['ligand_name']",
                        dest='input_path_dict',
                        type=str)
    parser.add_argument('--output_path',
                        help="Output path",
                        dest='output_path',
                        type=str)
    parser.add_argument('--dt',
                        help="Integration time in ps",
                        dest='dt',
                        default=0.004,
                        type=float)
    parser.add_argument('--simulation_time',
                        help="There are three different simulation: pulling, window equilibration, window production."
                        "You must provied the three corresponded times in this order (in ps). Defaults are:"
                        "pulling = 100000 ps"
                        "equilibration = 2000 ps"
                        "production = 20000 np",
                        dest='simulation_time',
                        nargs=3,
                        default=[100000, 2000, 20000],
                        type=float)
    parser.add_argument('--number_frame',
                        help="There are three different simulation: pulling, window equilibration, window production."
                        "You must provied the three corresponded number of frames in this order. Defaults are:"
                        "pulling = 2000"
                        "equilibration = 2"
                        "production = 2000",
                        dest='number_frame',
                        nargs=3,
                        default=[2000, 2, 2000],
                        type=int)
    parser.add_argument('--temperature',
                        help="Temperature in K",
                        dest='temperature',
                        default=303.15,
                        type=float)
    parser.add_argument('--pull_coord_vec',
                        help="This is the vector along the pulling will be done. Default is [0.0, 0.0, -1.0]",
                        dest='pull_coord_vec',
                        nargs=3,
                        default=[0.0, 0.0, -1.0],
                        type=float)
    parser.add_argument('--flat_bottom_init',
                        help="The radium for the flat bottom restrain. Default is 0.5 nm",
                        dest='flat_bottom_init',
                        default=0.5,
                        type=float)
    parser.add_argument('--flat_bottom_k',
                        help="The force constant for the flat bottom restrain. Default is 400 kJ/mol*nm^2",
                        dest='flat_bottom_k',
                        default=400,
                        type=float)
    parser.add_argument('--orient_restrain_init',
                        help="The angle for the orientation restrain. Default is 0.349 rad, 20 degrees 20/180*3.141",
                        dest='orient_restrain_init',
                        default=0.349,
                        type=float)
    parser.add_argument('--orient_restrain_k',
                        help="The force constant for the orientation restrain. Default is 1308.36 kJ/mol/rad^2, leads to sigma = 2.5 degrees",
                        dest='orient_restrain_k',
                        default=1308.36,
                        type=float)
    parser.add_argument('--orient_restrain_vec',
                        help="This means that the vector will be calculate in suing the initial orientation fo the ligands. if a vector is provied, for example (0,0,-1), then this will be used. You need to provied the vector in the format '0 0 -1'. 'calculate' is default",
                        dest='orient_restrain_vec',
                        default="calculate",
                        type=str)
    parser.add_argument('--pull_distance',
                        help="''",
                        dest='pull_distance',
                        default=4.5,
                        type=float)
    parser.add_argument('--pull_middle_points',
                        help="''",
                        nargs="+",
                        dest='pull_middle_points',
                        default=[],
                        type=float)
    parser.add_argument('--pull_force_constants',
                        help="''",
                        nargs="+",
                        dest='pull_force_constants',
                        default=[10000],
                        type=float)
    parser.add_argument('--window_widths',
                        help="''",
                        nargs="+",
                        dest='window_widths',
                        default=[0.1],
                        type=float)
    parser.add_argument('--chains',
                        help="''",
                        nargs="+",
                        dest='chains',
                        default=['A', 'B', 'C', 'D', 'E'],
                        type=str)
    parser.add_argument('--GROMACS_version',
                        help="",
                        dest='GROMACS_version',
                        default="2021.1",
                        type=str)
    parser.add_argument('--cpu',
                        help="",
                        dest='cpu',
                        default=12,
                        type=int)
    parser.add_argument('--COM_dist_cpu',
                        help="",
                        dest='COM_dist_cpu',
                        default=12,
                        type=int)
    parser.add_argument('--machine',
                        help='Type of machine smaug, gwdg or frontend',
                        dest='machine',
                        default="smaug",
                        type=str)
    parser.add_argument('--exclude_nodes',
                        help="",
                        dest='exclude_nodes',
                        default="",
                        type=str)
    parser.add_argument('--pull_code',
                        help='If you want to execute the pulling code add the flag . Default is False',
                        nargs="?",
                        dest='pull_code',
                        const=True,
                        default=False,
                        type=bool)
    parser.add_argument('--window_code',
                        help='If you want to execute the window code. Default is False',
                        nargs="?",
                        dest='window_code',
                        const=True,
                        default=False,
                        type=bool)
    args = parser.parse_args()

    args.input_path_dict = eval(args.input_path_dict)
    args.simulation_time = {
        'pulling': args.simulation_time[0], 'equilibration': args.simulation_time[1], 'production': args.simulation_time[2]}
    args.number_frame = {
        'pulling': args.number_frame[0], 'equilibration': args.number_frame[1], 'production': args.number_frame[2]}
    args.orient_restrain_vec = args.orient_restrain_vec.strip().lower()
    if args.orient_restrain_vec != "calculate":
        args.orient_restrain_vec = [
            float(xi) for xi in args.orient_restrain_vec.split()]

    main(args.input_path_dict,
         args.output_path,
         dt=args.dt,
         simulation_time=args.simulation_time,
         number_frame=args.number_frame,
         temperature=args.temperature,
         pull_coord_vec = args.pull_coord_vec,
         flat_bottom_init=args.flat_bottom_init,
         flat_bottom_k=args.flat_bottom_k,
         orient_restrain_init=args.orient_restrain_init,
         orient_restrain_k=args.orient_restrain_k,
         orient_restrain_vec=args.orient_restrain_vec,
         pull_distance=args.pull_distance,
         pull_middle_points=args.pull_middle_points,
         pull_force_constants=args.pull_force_constants,
         window_widths=args.window_widths,
         chains=args.chains,
         GROMACS_version=args.GROMACS_version,
         cpu=args.cpu,
         COM_dist_cpu=args.COM_dist_cpu,
         machine=args.machine,
         exclude_nodes=args.exclude_nodes,
         pull_code=args.pull_code,
         window_code=args.window_code)
