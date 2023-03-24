#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
The correct way to use this script is from the above directory and
htmd, to activate the environment
nohup ./auto_umbrella.py >output.txt 2>&1 &
"""
from palmiche.umbrella import make_posre, COM_distance, conf_windows
from palmiche.utils import tools, mdp, jobsh, xvg, markov
from palmiche.simultemp import weights
import os, datetime, tempfile, argparse, time, multiprocessing, socket
import numpy as np
import pandas as pd

def dict_gen(path = "./HMR_2_5", FF_protein = "amber99sb-star-ildn.ff",
             FF_lipid="Slipids_2020.ff", FF_ligand="GAFF2.ff", topol="topol_heavyH.top",
             toppar="toppar", conf="confout.gro"):
    """
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
                                                           "conf": os.path.join(path, receptor, ligand, conf),
                                                           }
    return dict_paths

def setup_production(
    frame,
    frame2window,
    log_center_cluster_abspath,
    simulated_tempering_weights_equilibration_time,
    new_path_dict,
    mdp_options,
    simulation_time,
    number_frame,
    pull_distance,
    ligands,
    flat_bottom_init,
    flat_bottom_k,
    orient_restrain,
    orient_restrain_init,
    orient_restrain_k,
    orient_restrain_vec,
    dt,
    windows,
    pull_coord_vec,
    annealing,
    tc_grps,
    temperature,
    exclude_nodes,
    cpu,
    partition,
    GROMACS_version,
    str_simulated_tempering_alpha,
    simulated_tempering_temperatures,
    frontend,
    mdrun_opt,
    ):
    # Esto tegno que ponerlo dentro de un civlo while qu espere y que duerma por un numero rando de tiempo entre 10 y 120, para que se puedan desfasar los threads,
    #jobs2launch[ID_dir_abspath] = {'cluster': cluster, 'window_ID':window_ID}
    while True:
        # Check if the log file exist
        # Esto ciclo es peligroso, porque en caso de errores el caclulo se quedara enganchado esperadno a que aparezca el production.log file
        # Aqui seria bueno algo que me de un terminado en caso de que pase mucho y no suceda nada talvel puedo medir el tiempo y si es mayor que eso y cae en esa ocndicion rompo
        if os.path.isfile(log_center_cluster_abspath):
            time_progression, weights_info = weights.get_weights_from_log(log_center_cluster_abspath, plot = False)
            # Check if the time for the equilibration of the weights was reached
            if time_progression[-1] >= simulated_tempering_weights_equilibration_time:
                # Prepearing all the necesary files and directories for the memebers of the current cluster

                window_ID = frame2window[frame]
                conf = f"conf{frame}.gro"
                ID_dir_basename = str(window_ID).zfill(5)
                tools.makedirs(ID_dir_basename)
                ID_dir_abspath_member = os.path.abspath(ID_dir_basename)
                print(f'I am inside of setup_production {ID_dir_abspath_member}')
                tools.mv(os.path.join('split_xtc',conf), ID_dir_abspath_member)
                [tools.cp(new_path_dict[key], ID_dir_abspath_member, r=True) for key in new_path_dict]
                #   Create mdp, equilibration
 
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
                    orient_restrain=orient_restrain,
                    orient_restrain_init=orient_restrain_init,
                    orient_restrain_k=orient_restrain_k,
                    orient_restrain_vec=orient_restrain_vec,
                    # General options for the actual pulling
                    pull_nstxout = int(simulation_time['equilibration'] / (dt * number_frame['equilibration'])), # Reduce the output frequency of nstxout. It is not relevant. The frequency was set equal to the nstxout_compressed
                    pull_nstfout = int(simulation_time['equilibration'] / (dt * number_frame['equilibration'])), # Reduce the output frequency of nstfout. It is not relevant. The frequency was set equal to the nstxout_compressed
                    pull_coord1_k=windows[window_ID]['force_constant'],
                    pull_coord1_vec=" ".join([str(xi) for xi in pull_coord_vec]),
                    pull_coord1_init=windows[window_ID]['dist'],
                    pull_coord1_rate= 0.0,)

                # Adding the section for the annealing if required.
                if annealing:
                    # This is a fix method
                    # Increase in 50 K the temperature for the ligands and the protein
                    # during the initial 25 % of the time, then the continued with a cooling (50%)
                    # And finally set constant the temperature to the initial one
                    GroupBool = (np.array(tc_grps) != 'MEMB') & (np.array(tc_grps) != 'SOLV')
                    temp_list = np.arange(temperature, temperature + 52, 2)
                    MDP_equilibration.annealing2(GroupBool, temp_list, heating_frac=0.25, cooling_frac=0.50)
                MDP_equilibration.write(os.path.join(ID_dir_abspath_member, 'equilibration.mdp'))

                #   Create mdp production with the corresponded weights
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
                    orient_restrain=orient_restrain,
                    orient_restrain_init=orient_restrain_init,
                    orient_restrain_k=orient_restrain_k,
                    orient_restrain_vec=orient_restrain_vec,
                    # General options for the actual pulling
                    pull_coord1_k=windows[window_ID]['force_constant'],
                    pull_coord1_vec=" ".join([str(xi) for xi in pull_coord_vec]),
                    pull_coord1_init=windows[window_ID]['dist'],
                    pull_coord1_rate= 0.0,)
                MDP_production.simulated_tempering(
                    temp_list = simulated_tempering_temperatures, # The same used for the anealing for the prediction of the inital weights
                    init_lambda_weights = weights_info[-1,:,1], # Take final weights and assign them to their respective members of cluster
                )
                MDP_production.write(os.path.join(ID_dir_abspath_member,'production.mdp'))
                # Create the 3_job.sh, for the simulated tmeepring steps (production) -update gpu can not be used becasue is incompatible with the integration method md-vv (the one used for ST)
                JOBSH_production = jobsh.JOB(
                    sbatch_keywords={
                        'job-name':f"{window_ID}_production",
                        'exclude':exclude_nodes,
                        'cpus-per-task': cpu,
                        'partition': partition,
                    },
                    GROMACS_version = GROMACS_version,
                    build_GROMACS_section=f"{str_simulated_tempering_alpha}\n"\
                        f"gmx grompp -f equilibration.mdp -o equilibration.tpr -c {conf} -r {conf} -p {os.path.basename(new_path_dict['topol'])} -n {os.path.basename(new_path_dict['ndx'])}\n"\
                        f"gmx mdrun -nt {cpu} -cpi -stepout 5000 -v -deffnm equilibration -px equilibration_pullx -pf equilibration_pullf {mdrun_opt} >& equilibration.lis\n\n\n"\
                        f"gmx grompp -f production.mdp -o production.tpr -c equilibration.gro -r equilibration.gro -p {os.path.basename(new_path_dict['topol'])} -n {os.path.basename(new_path_dict['ndx'])}\n"\
                        f"gmx mdrun -nt {cpu} -cpi -stepout 5000 -v -deffnm production -px production_pullx -pf production_pullf {mdrun_opt} >& production.lis\n"
                )
                JOBSH_production.write(os.path.join(ID_dir_abspath_member,'3_job.sh'))

                # Launch the jobs and wait for end of execution
                print(f'we are launching inside of the setup_production() the production job for the corresponded cluster memebers. Jobs = {ID_dir_abspath_member}')
                if frontend:
                    tools.job_launch_list([os.path.join(ID_dir_abspath_member, '3_job.sh')], shell='bash')# This is just for testing
                else:
                    # This will launch and reapeat (jobsh_name) in case that logfile did not finish correctlly due to runing out of tiem on the selected partition,
                    tools.launch_wait_check_repeat(
                        partition = partition,
                        jobpaths = [ID_dir_abspath_member],
                        logfile = 'production.log',
                        lisfile = 'production.lis',
                        jobsh_name = '3_job.sh',
                        logfile_start_datetime = 'equilibration.log',)
                # If the condition was reached and all the jobs were launched and correctlly checked, breack the infinite while loop
                print(f"The launching inside of the setup_production() the production of the Jobs = {os.path.join(ID_dir_abspath_member, '3_job.sh')} was succesed.")
                break
            # If the time has not been reached.  sleep random amount of seconds (this randomness is used to defases the threads)
            else:
                time.sleep(np.random.randint(20,60))
        # If log file doesn't exist.  sleep random amount of seconds (this randomness is used to defases the threads)
        else:
            time.sleep(np.random.randint(20,60))

def main(input_path_dict,
        output_path,
        dt=0.004,
        simulation_time={'pulling': 100000,
                        'equilibration': 2000,
                        'annealing_for_ST':5000, # For now this is irrelevant because we are using the PandeWeights and for that the methodology is already fixed. This will be usefull in the case of the use of AleWeights
                        'production': 20000}, # ps
        number_frame={'pulling': 2000,
                    'equilibration': 2,
                    'annealing_for_ST': 20,
                    'production': 2000},
        temperature = 303.15,
        refine_init_pull = False,
        annealing = True,
        pull_coord_vec = [0.0, 0.0, -1.0],
        reverse = False, # If True, it will use the negative pull_coord1_rate, Will go in the reverse direction with the negative of the reaction coordinates.
        flat_bottom_init=0.4,  # This is the radium
        flat_bottom_k=500, # This leads to sigma = 0.07942668639322366, this give us nice results, but we could test with 1000 that lead around 0.05 sigma
        orient_restrain = True,
        lig_vec_group1 = None,
        lig_vec_group2 = None,
        orient_restrain_init=45,  # 0.349, #; 20 degrees 20/180*3.141
        orient_restrain_k=2.5,  # #; in kJ/mol/rad^2, leads to sigma = 1 degrees k = KbT/sigma^2
        # This means that the vector will be calculate in suing the initial orientation fo the ligands. if a vector is provied, for example (0,0,-1), then this will be used
        orient_restrain_vec="calculate",
        pull_distance=4.5,
        # Here I will identify what are the different intervals on the windows.
        pull_middle_points=None,
        # Here I will deffine how many pulling force constant I will use, if there is only on number then the same applied for every thong, I could even add the pulling force constant for the pulling itself and the one for the windows
        pull_force_constants=[10000],
        # This variable also will have as many as different intervals are specified, in case only one, tha same applied for all the windows
        window_widths=[0.1],
        chains=['A', 'B', 'C', 'D', 'E'],
        # The portein is always the first group and has as name SOLU, so it is fine to create a selection based only on the protein and
        # should work for different ligands as well.
        CLOSE_AA_radius = 0.4, # This is for the selction of the close aa, by default is 0.4 nm 
        CLOSE_AA_NDX = None, # This is a GROMACS ndx file it must have as many group as chains and the group names must be LI{chain}_CLOSE_AA
        GROMACS_version="2021.4",
        cpu=12,
        COM_dist_cpu=12,
        frontend=False,
        partition="deflt",
        exclude_nodes="",
        pull_code=True,
        window_code=True,
        # Here I will define the variables for the simulated tampering in the windows
        # If True weights.PandeWeights will be the method to use, in the future we should add our method and compare.
        # The only parameter that we will controll is alpha
        # The spacing between temperatures will be fixed to 2
        # And the maximum temperature will be fixed to 30K higher that the value of temperatures
        # The cluster method will be fpi (we could also use the cheaper one distance and number of interactions)
        # If annealing is True, first an equilibration will be done with the already implemented annealing and after that
        # A second annealing will be done for getting the total_energies over the last structure obtained for the previus annealing
        # You need to provided the finger print file with name fp_df.pkl
        simulated_tempering = True,
        simulated_tempering_temperatures = np.arange(303.15, 303.15 + 32, 2), # A sorted list of temperatures, the first temperature must be the temperature of provided previously.
        simulated_tempering_alpha = None,#0.5, Only useful if 2021.6_simTemp_setOccupState0 GROMACS'version is used (only available on smaug at the moment)
        simulated_tempering_n_clusters = 10,
        simulated_tempering_weights_equilibration_time = 10000, #in ps. For now this option is only compatible with simulated_tempering = True and simulated_tempering_n_clusters = 0
        test = False,
        ):
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
        DESCRIPTION. The default is {'pulling':100000,'equilibration':1000,'production':10000}. [ps]
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

    Tengo que completar la documentacion
    Returns
    -------
    None.

    """
    # This check is for the case of testing simulations
    if frontend:
        hostname = socket.gethostname()
        update_gpu = ''
        if hostname == 'smaug':
            mdrun_opt = '-pin on -pinstride 1 -pinoffset 12 -gpu_id 1 -nice 19'# '-pin on -pinstride 1 -pinoffset 0 -gpu_id 0 -nice 19'
        else:
            mdrun_opt = '-pin on -pinstride 1 -pinoffset 0 -nice 19'
    else:
        mdrun_opt = ''
        update_gpu = '-update gpu'

    # Check for the ST temperatures list
    if simulated_tempering:
        simulated_tempering_temperatures = np.sort(simulated_tempering_temperatures)
        if simulated_tempering_temperatures[0] != temperature:
            raise ValueError(f"The lower temperature of simulated_tempering_temperatures = {simulated_tempering_temperatures[0]} must be the same as temperature = {temperature}")

    # See if the alpha was selected for ST and creating the strin for the jobsh
    if simulated_tempering_alpha:
       str_simulated_tempering_alpha = f"export SIMTEMP_OCCUP_STATE0={simulated_tempering_alpha}\n"
       if GROMACS_version != '2021.6_simTemp_setOccupState0':
           raise ValueError(f"simulated_tempering_alpha = {simulated_tempering_alpha} but GROMACS_version = {GROMACS_version}. Only GROMACS_version = 2021.6_simTemp_setOccupState0 is possible in this case")
    else:
        str_simulated_tempering_alpha = '\n'
    # Here I am see how many different intervals are in the windows
    # Selecting the correct set of parameters for pull_force_constants and window_widths
    # This is important, it will always sorted pull_midles_points, and pull_force_constants and window_widths it will be map respect to the sorted pull_middle_points
    # I have to consider the case when non pull_middle_points are set. In this case I also have to chane some other option when the windows are building.
    pull_middle_points = sorted(pull_middle_points)

    # Check that pull_coord_vec is correect
    if pull_coord_vec[0]!=0 or pull_coord_vec[1] !=0:
        raise ValueError(f"pull_coord_vec must be a vector (0,0,something) not {pull_coord_vec}")
    
    # Generating the ligands name whit the chain information
    ligands = [f"LI{chain}" for chain in chains]
    # Create grouplist in the correct order
    grouplist = [[ligand+"_CLOSE_AA",ligand] for ligand in ligands]
    # # Swap the order if needed to keep the correct sign on the distnce
    # if (reverse and pull_coord_vec[-1] > 0) or (not reverse and pull_coord_vec[-1] < 0):
    #     for group in grouplist:
    #         group[0], group[1] = group[1], group[0]
    
    print(
        f"The following window intervals will be constructed:\nmin---{'---'.join([f'[{point}]' for point in pull_middle_points])}---max\n")
    if len(pull_force_constants) == 1:
        pull_force_constants = (
            len(pull_middle_points) + 2) * pull_force_constants
        print(
            f"The force constant {pull_force_constants[0]} will be used for all the simulations.\n")
    elif len(pull_force_constants) == 2:
        pull_force_constants = [pull_force_constants[0]] + [pull_force_constants[1] for i in range(len(pull_middle_points) + 1)]
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

    cwd = os.getcwd()
    start_datetime = datetime.datetime.now()
    output_path = os.path.abspath(output_path)

    # Generating the tc_grps names
    tc_grps = ['SOLU', 'MEMB', 'SOLV'] + ligands
    # The general options for the mdp
    mdp_options = {
        "dt": dt,
        "tc_grps": ' '.join(tc_grps),
        "tau_t": ' '.join(len(tc_grps)*['1.0']),
        "ref_t": ' '.join(len(tc_grps)*[str(temperature)])
    }
    # Here I deffine the variable names for the index.
    out_index_name = 'index.ndx'

    window_path = os.path.join(output_path, 'windows')
    new_path_dict = {}
    for key in input_path_dict:
        if key != 'conf': # I don't want to copy this file to the directory, so, I need to keep the original path.
            new_path_dict[key] = os.path.join(
                output_path, os.path.basename(input_path_dict[key]))

    new_path_dict['ndx'] = os.path.join(output_path, 'index.ndx')
    [tools.makedirs(item) for item in [output_path, window_path]]

    # The two if for pull_code and for window_code is in order to flexibly the running of the script
    if pull_code:
        # Here I am creating the groups that I need, The ones for the pulling code.
        # And the rest are for the temperature coupling.

        # Creating the tmp files
        tmpopt1 = tempfile.NamedTemporaryFile(suffix='.opt')
        tmpopt2 = tempfile.NamedTemporaryFile(suffix='.opt')
        tmpndx1 = tempfile.NamedTemporaryFile(suffix=".ndx")
        tmpndx2 = tempfile.NamedTemporaryFile(suffix=".ndx")
        tmpndx3 = tempfile.NamedTemporaryFile(suffix=".ndx")
        tmpndx4 = tempfile.NamedTemporaryFile(suffix=".ndx")

        tools.run(
            f"export GMX_MAXBACKUP=-1; echo q | gmx make_ndx -f {input_path_dict['conf']} -o {tmpndx1.name}")
        options = f"""
        "SOLU" group Protein;
        "MEMB" ((group System and ! group Water_and_ions) and ! group Protein) and ! (resname {" or resname ".join(ligands)});
        "SOLV" group Water_and_ions;
        "Backbone" group Backbone;
        """
        if CLOSE_AA_NDX:
            if len(tools.NDX(CLOSE_AA_NDX).data) != len(chains):
                CLOSE_AA_NDX = None
                print(f"CLOSE_AA_NDX have {len(CLOSE_AA_NDX)} index groups, but you have {len(chains)} chains. The index groups will be automatically genereated.")

        for ligand in ligands:
            options += f"\"{ligand}\" resname {ligand};\n"
            if not CLOSE_AA_NDX:
                options += f"\"{ligand}_CLOSE_AA\" group \"Backbone\" and same residue as within {CLOSE_AA_radius} of resname {ligand};\n"

        # ==============================================================================

        with open(tmpopt1.name, 'w') as f:
            f.write(options)
        tools.run(
            f"export GMX_MAXBACKUP=-1; gmx select -s {input_path_dict['conf']} -sf {tmpopt1.name} -n {tmpndx1.name} -on {tmpndx2.name}")
        # deleting the the line _f0_t0.000 in the file
        with open(tmpndx2.name, 'r') as f:
            data = f.read()
            data = data.replace("_f0_t0.000", "")
        with open(tmpndx2.name, "w") as f:
            f.write(data)

        # Rectification, to take the same aa in all the monomers
        if CLOSE_AA_NDX:
            close_ndx = tools.NDX(CLOSE_AA_NDX).data
            ndx_rect_dict = tools.NDX(tmpndx2.name).data
            ndx_rect_dict.update(close_ndx)
        else:
            ndx_rect_dict = make_posre.ndx_rectification(tmpndx2.name, input_path_dict['conf'], [
                                                    f"{ligand}_CLOSE_AA" for ligand in ligands])

        # Adding the corresponding index groups for the orientation restrain.
        if orient_restrain:
            for ligand in ligands:
                ndx_rect_dict[f"{ligand}_GROUP1"] = [ndx_rect_dict[ligand][index - 1]
                                                for index in lig_vec_group1]
                ndx_rect_dict[f"{ligand}_GROUP2"] = [ndx_rect_dict[ligand][index - 1]
                                                for index in lig_vec_group2]
        ndx_obj = tools.NDX(tmpndx3.name)
        ndx_obj.data = ndx_rect_dict
        ndx_obj.write(out_name=tmpndx3.name,backup=False)

        # to rectify groups, just keep the atoms that belong to the group backbone
        # , the rest is just repeat to take all the groups
        options = """
        "SOLU" group SOLU;
        "MEMB" group MEMB;
        "SOLV" group SOLV;
        """
        for i, ligand in enumerate(ligands):
            options += f"\"{ligand}\" group {ligand};\n"\
            f"\"{ligand}_CLOSE_AA\" group {ligand}_CLOSE_AA and group \"Backbone\";\n"

            if orient_restrain:
                options += f"\"{ligand}_GROUP1\" group {ligand}_GROUP1;\n"\
                f"\"{ligand}_GROUP2\" group {ligand}_GROUP2;\n"

        with open(tmpopt2.name, 'w') as f:
            f.write(options)
        tools.run(
            f"export GMX_MAXBACKUP=-1; gmx select -s {input_path_dict['conf']} -sf {tmpopt2.name} -n {tmpndx3.name} -on {tmpndx4.name}")
        with open(tmpndx4.name, 'r') as f:
            data = f.read()
            data = data.replace("_f0_t0.000", "")
        with open(tmpndx4.name, "w") as f:
            f.write(data)

        # Ordering the files, now all the path in new_path_dict are useful
        [tools.cp(input_path_dict[key], output_path, r=True) for key in input_path_dict if key!= 'conf']
        tools.cp(tmpndx4.name, os.path.join(output_path, out_index_name))

        # Changing the directory to
        os.chdir(output_path)
        # Creating the average oriented vector. This vector is used for the orientation restrain of the ligands. This is done is order to have
        # more consisten the simulation and the stadistic better for the 5 PMF, because if not he ligands have different orientation, and that means that the PMFs are not comparables.
        # I will calculate only fi the vector provi
        if orient_restrain:
            try:  # I need to use try, because if a list was provided then, will raise an error, because the list doesn't have strip() method
                if orient_restrain_vec.strip().lower() == "calculate":
                    # first I need to create a temporal tpr file
                    mdptmp = tempfile.NamedTemporaryFile(suffix=".mdp")
                    tprtmp = tempfile.NamedTemporaryFile(suffix=".tpr")
                    tools.run(
                        f"export GMX_MAXBACKUP=-1; gmx grompp -f {mdptmp.name} -c {os.path.basename(input_path_dict['conf'])} -r {os.path.basename(input_path_dict['conf'])} -p {os.path.basename(input_path_dict['topol'])} -n {out_index_name} -o {tprtmp.name}")

                    # Here I am calculating all the vectors for each ligand
                    vectors = [tools.get_vec_COM(input_path_dict['conf'], tprtmp.name,
                                                out_index_name, f"{ligand}_GROUP1", f"{ligand}_GROUP2") for ligand in ligands]
                    # Get the average oriented vector, rounded with 3 decimals
                    orient_restrain_vec = tools.aovec(vectors, round=3)
            except:
                pass
        # Generating the mdp using the templates module, here we are constructing the pulling.
        pull_coord1_init_pulling = 0
        if refine_init_pull:
            # Getting the initial distnace to avoid future clashes
            print('Building the refinment for the initial configuaration of the pulling')
            mdp_pull_tmp = tempfile.NamedTemporaryFile(suffix=".mdp")
            tpr_pull_tmp = tempfile.NamedTemporaryFile(suffix=".tpr")
            tools.run(f"export GMX_MAXBACKUP=-1; gmx grompp -f {mdp_pull_tmp.name} -c {input_path_dict['conf']} -r {input_path_dict['conf']} -p {input_path_dict['topol']} -n {out_index_name} -o {tpr_pull_tmp.name}")
            pull_coord1_init_pulling = COM_distance.COM_dist_chunk_code(
                conf_file = input_path_dict['conf'],
                grouplist = grouplist,
                ndx =out_index_name,
                tpr = tpr_pull_tmp.name,
                prefix = None,
                projected_on = 'z')[1]

            # Creating the pull_minimization mdp
            MDP_pulling_minimization = mdp.MDP(
                type = 'minimization',
                define = f"-DPOSRES -DPOSRES_FC_BB=1000.0 -DPOSRES_FC_SC=100.0 -DPOSRES_FC_LIPID=1000.0 -DDIHRES -DDIHRES_FC=1000.0 -DPOSRES_LIG=0.0"
                )
            MDP_pulling_minimization.pull(
                pull_distance,
                ligands,
                # Options for the cylinder and orientation restraints
                flat_bottom_init=flat_bottom_init,
                flat_bottom_k=flat_bottom_k,
                orient_restrain=orient_restrain,
                orient_restrain_init=orient_restrain_init,
                orient_restrain_k=orient_restrain_k,
                orient_restrain_vec=orient_restrain_vec,
                # General options for the actual pulling
                pull_nstxout = 0,
                pull_nstfout = 0,
                pull_coord1_rate = 0.0,
                pull_coord1_k=pull_force_constants[0],
                pull_coord1_vec=" ".join([str(xi) for xi in pull_coord_vec]),
                pull_coord1_init=pull_coord1_init_pulling,
            )
            MDP_pulling_minimization.write(os.path.join(output_path, 'pull_minimization.mdp'))

            # Creating the pull_equilibration mdp
            # Reduce the time steep
            mdp_options_copy  = mdp_options.copy()
            mdp_options_copy['dt'] = 0.002
            if test:
                time_pulling_equilibration = 50
            else:
                time_pulling_equilibration = 1000
            MDP_pulling_equilibration = mdp.MDP(
                    **mdp_options_copy,
                    define = f"-DPOSRES -DPOSRES_FC_BB=1000.0 -DPOSRES_FC_SC=100.0 -DPOSRES_FC_LIPID=1000.0 -DDIHRES -DDIHRES_FC=1000.0 -DPOSRES_LIG=0.0",
                    time = time_pulling_equilibration,
                    xtc_numb_frame = 5
            )
            MDP_pulling_equilibration.pull(
                pull_distance,
                ligands,
                # Options for the cylinder and orientation restraints
                flat_bottom_init=flat_bottom_init,
                flat_bottom_k=flat_bottom_k,
                orient_restrain=orient_restrain,
                orient_restrain_init=orient_restrain_init,
                orient_restrain_k=orient_restrain_k,
                orient_restrain_vec=orient_restrain_vec,
                # General options for the actual pulling
                pull_nstxout = 0,
                pull_nstfout = 0,
                pull_coord1_rate = 0.0,
                pull_coord1_k=pull_force_constants[0],
                pull_coord1_vec=" ".join([str(xi) for xi in pull_coord_vec]),
                pull_coord1_init=pull_coord1_init_pulling,
            )
            MDP_pulling_equilibration.write(os.path.join(output_path, 'pull_equilibration.mdp'))

        print('Building the pulling mdp')
        MDP_pulling = mdp.MDP(
            **mdp_options,
            define = f"-DPOSRES -DPOSRES_FC_BB=1000.0 -DPOSRES_FC_SC=100.0 -DPOSRES_FC_LIPID=1000.0 -DDIHRES -DDIHRES_FC=1000.0 -DPOSRES_LIG=0.0",
            time = simulation_time['pulling'],
            xtc_numb_frame = number_frame['pulling'])
        pull_coord1_rate = pull_distance / simulation_time['pulling']
        if reverse: pull_coord1_rate = -pull_coord1_rate
        MDP_pulling.pull(
            pull_distance,
            ligands,
            # Options for the cylinder and orientation restraints
            flat_bottom_init=flat_bottom_init,
            flat_bottom_k=flat_bottom_k,
            orient_restrain=orient_restrain,
            orient_restrain_init=orient_restrain_init,
            orient_restrain_k=orient_restrain_k,
            orient_restrain_vec=orient_restrain_vec,
            # General options for the actual pulling
            pull_nstxout = int(simulation_time['pulling'] / (dt * number_frame['pulling'])), # Reduce the output frequency of nstxout. It is not relevant. The frequency was set equal to the nstxout_compressed
            pull_nstfout = int(simulation_time['pulling'] / (dt * number_frame['pulling'])), # Reduce the output frequency of nstfout. It is not relevant. The frequency was set equal to the nstxout_compressed
            pull_coord1_k=pull_force_constants[0],
            pull_coord1_vec=" ".join([str(xi) for xi in pull_coord_vec]),
            pull_coord1_init=pull_coord1_init_pulling,
            pull_coord1_rate = pull_coord1_rate,
        )
        MDP_pulling.write(os.path.join(output_path, 'pull.mdp'))

        if refine_init_pull:
            build_GROMACS_section_pull=f"{str_simulated_tempering_alpha}\n"\
                f"gmx grompp -f pull_minimization.mdp -c {input_path_dict['conf']} -r {input_path_dict['conf']} -p {os.path.basename(input_path_dict['topol'])} -n {out_index_name} -o pull_minimization.tpr\n"\
                f"gmx mdrun -nt {cpu} -cpi -stepout 5000 -v -deffnm pull_minimization -px pull_minimization_pullx -pf pull_minimization_pullf {mdrun_opt} >& pull_minimization.lis\n"\
                f"gmx grompp -f pull_equilibration.mdp -c pull_minimization.gro -r pull_minimization.gro -p {os.path.basename(input_path_dict['topol'])} -n {out_index_name} -o pull_equilibration.tpr\n"\
                f"gmx mdrun -nt {cpu} -cpi -stepout 5000 -v -deffnm pull_equilibration -px pull_equilibration_pullx -pf pull_equilibration_pullf {mdrun_opt} {update_gpu} >& pull_equilibration.lis\n"\
                f"gmx grompp -f pull.mdp -c pull_equilibration.gro -r pull_equilibration.gro -p {os.path.basename(input_path_dict['topol'])} -n {out_index_name} -o pull.tpr\n"\
                f"gmx mdrun -nt {cpu} -cpi -stepout 5000 -v -deffnm pull -px pull_pullx -pf pull_pullf {mdrun_opt} {update_gpu} >& pull.lis\n"
            logfile_start_datetime_pull = 'pull_minimization.log'
        else:
            build_GROMACS_section_pull=f"{str_simulated_tempering_alpha}\n"\
                f"gmx grompp -f pull.mdp -c {input_path_dict['conf']} -r {input_path_dict['conf']} -p {os.path.basename(input_path_dict['topol'])} -n {out_index_name} -o pull.tpr\n"\
                f"gmx mdrun -nt {cpu} -cpi -stepout 5000 -v -deffnm pull -px pull_pullx -pf pull_pullf {mdrun_opt} {update_gpu} >& pull.lis\n"
            logfile_start_datetime_pull = 'pull.log'

        JOBSH_pull = jobsh.JOB(
            sbatch_keywords={
                'job-name':'pull',
                'exclude':exclude_nodes,
                'cpus-per-task': cpu,
                'partition': partition,
            },
            GROMACS_version = GROMACS_version,
            build_GROMACS_section=build_GROMACS_section_pull,
        )
        JOBSH_pull.write()
        print('Launching pulling simulation')
        # Executing the job
        if frontend:
            tools.run(f"bash job.sh")
        else:
            # This will launch and reapeat (jobsh_name) in case that logfile did not finish correctlly due to runing out of tiem on the selected partition,
            tools.launch_wait_check_repeat(
                partition = partition,
                jobpaths = ['.'],
                logfile = 'pull.log',
                lisfile = 'pull.lis',
                jobsh_name = 'job.sh',
                logfile_start_datetime=logfile_start_datetime_pull)
            check = tools.CHECK("pull.log", "pull.lis")
            if not check.performance:
                raise RuntimeError(f"The pull simulation ended with error. Chek the simulation on {os.getcwd()}")
        print('The pulling simulation finished')

    if window_code:
        print('Starting the configuration of the windows')
        # Here We calculate the fingerprint if simulated tempering is turn on and simulated_tempering_n_clusters is set to a number different to zero.
        if simulated_tempering and simulated_tempering_n_clusters:
            os.chdir(output_path)
            protein_selection = 'segid *Protein_chain_A'
            ligand_selection = 'resname LIA'
            if os.path.isfile('ifp.pkl'): # Check if the interaction fingerprint exist
                classifier_ifp = tools.classifier_ifp('pull.tpr', 'pull.xtc',ligand_selection=ligand_selection, protein_selection=protein_selection, n_clusters =  simulated_tempering_n_clusters, load_fingerprint = True, ifp_file = 'ifp.pkl', ifp_atoms_file = 'ifp_atoms.pkl')
            else:
                classifier_ifp = tools.classifier_ifp('pull.tpr', 'pull.xtc', ligand_selection=ligand_selection, protein_selection=protein_selection, n_clusters =  simulated_tempering_n_clusters, load_fingerprint = False, ifp_file = 'ifp.pkl', ifp_atoms_file = 'ifp_atoms.pkl')

        # Changing dir to windows and measurement of distances
        os.chdir(window_path)

        # This will output a file called mean.dat, that has the format needed for pulling_processing.windows_creator. The distance is the average distance of the grouplist
        tools.makedirs("split_xtc")
        print('Calculating distances')
        COM_distance.main(grouplist, cpu=COM_dist_cpu, ndx=new_path_dict['ndx'], tpr=os.path.join(output_path,'pull.tpr'),
                        xtc=os.path.join(output_path, 'pull.xtc'), prefix='conf', out='summary_distances.dat',
                        # Doing check of the distance "sign" and specifiying the columns
                        projected_on = 'z',
                        split_out_dir='split_xtc')
        # ==============================================================================
        # The processing and the launch of the windows

        # ==============================================================================

        # _______checking__________
        #check_cont = input("Now the configuration of the windows will take place. Will you continue?[y/n]:").lower()[0]
        #check_cont = "y"
        # _________________________
        # if check_cont == "y":
        # This part of the code is for the creation of the windows with their respective pulling forces and interval widths.
        print('Selection the frames for each window')
        windows = {}
        # I should set minimum=pull_middle_points[i - 1] + width, but think on it first
        if pull_middle_points:
            maximum_windows_ID = 0
            for (i, width) in enumerate(window_widths):
                if i == 0:
                    tmp_windows = conf_windows.main("mean.dat", width, txt_summary=None, maximum=pull_middle_points[i])
                    for window_ID in tmp_windows:
                        windows[window_ID] = tmp_windows[window_ID]
                        windows[window_ID]['interval'] = f"min-{pull_middle_points[i]}"
                        windows[window_ID]['width'] = width
                        # Because the first force_constant is for the pulling
                        windows[window_ID]['force_constant'] = pull_force_constants[i+1]
                elif i == len(window_widths) - 1:
                    # Here is added by default the sampled distance to the maximum value
                    tmp_windows = conf_windows.main(
                        "mean.dat", width, txt_summary=None,  minimum=pull_middle_points[i - 1])

                    for window_ID in tmp_windows:
                        # Only if the frame is not already in the dict, I will add the window
                        if tmp_windows[window_ID]['frame'] not in [windows[window_ID2]['frame'] for window_ID2 in windows]:
                            windows[window_ID + maximum_windows_ID] = tmp_windows[window_ID]
                            windows[window_ID + maximum_windows_ID]['interval'] = f"{pull_middle_points[i - 1]}-max"
                            windows[window_ID + maximum_windows_ID]['width'] = width
                            # Because the first force_constant is for the pulling
                            windows[window_ID + maximum_windows_ID]['force_constant'] = pull_force_constants[i+1]
                        else:
                            # If the frame is repeated I need to substract the len of the tmp_window
                            maximum_windows_ID -= 1
                else:
                    tmp_windows = conf_windows.main("mean.dat", width, txt_summary=None,  minimum=pull_middle_points[i - 1], maximum=pull_middle_points[i])

                    for window_ID in tmp_windows:
                        # Only if the frame is not already in the dict, I will add the window
                        if tmp_windows[window_ID]['frame'] not in [windows[window_ID2]['frame'] for window_ID2 in windows]:
                            windows[window_ID + maximum_windows_ID] = tmp_windows[window_ID]
                            windows[window_ID +
                                    maximum_windows_ID]['interval'] = f"{pull_middle_points[i - 1]}-{pull_middle_points[i]}"
                            windows[window_ID + maximum_windows_ID]['width'] = width
                            # Because the first force_constant is for the pulling
                            windows[window_ID + maximum_windows_ID]['force_constant'] = pull_force_constants[i+1]
                        else:
                            # If the frame is repeated I need to substract the len of the tmp_window
                            maximum_windows_ID -= 1
                maximum_windows_ID += len(tmp_windows)
        else:
            windows = conf_windows.main("mean.dat", window_widths[0], txt_summary=None)
            for window_ID in windows:
                windows[window_ID]['interval'] = "min-max"
                windows[window_ID]['width'] = window_widths[0]
                # Because the first force_constant is for the pulling
                windows[window_ID]['force_constant'] = pull_force_constants[1]

        # Creating converting dictionary frame to window window_ID
        frame2window = {windows[window_ID]['frame']:window_ID for window_ID in windows}

        if simulated_tempering:
            print('Configurating  ST. Cluster center and memebership')
            if simulated_tempering_n_clusters:
                # Here We get the membership of each window
                user_membership = {}
                for window_ID in windows:
                    windows[window_ID]['cluster'] = classifier_ifp.model.labels_[windows[window_ID]['frame']]
                    if windows[window_ID]['cluster'] in user_membership:
                        user_membership[windows[window_ID]['cluster']].append(windows[window_ID]['frame'])
                    else:
                        user_membership[windows[window_ID]['cluster']] = [windows[window_ID]['frame']]

                frame_close_to_centroids = classifier_ifp.get_closers(user_membership=user_membership)
                # We eliminate from user_membership the centers of clusters, this will be used apart.
                # frame_close_to_centroids is a dictionary with keys the cluster identifier and values the closest frame to the respective cluster centroid.
                for cluster in user_membership:
                    user_membership[cluster].remove(frame_close_to_centroids[cluster])
                # Now we have:
                #   frame_close_to_centroids; that will be used to know what windows will be used for the initial guessing of weights
                #   user_membership; to assign the corresponded initial weights from the final simulation of the windows belongs to frame_close_to_centroids
            else:
                # All will be in the same cluster (0)
                # The window used for the calculation will be window 1
                cluster = 0
                window_ID_centroid = 1
                user_membership = {cluster:[]}
                for window_ID in windows:
                    windows[window_ID]['cluster'] = 0
                    user_membership[cluster].append(windows[window_ID]['frame'])

                frame_close_to_centroids = {cluster: windows[window_ID_centroid]['frame']}
                user_membership[cluster].remove(frame_close_to_centroids[cluster])

        print('Saving window configuration data')
        # Construction of the DataFrame to print
        df = pd.DataFrame.from_dict(windows, orient='index')
        df.insert(0, 'window_ID', df.index)
        df.reset_index(drop=True, inplace=True)
        # Rectifying the delta_dit column
        df.delta_dist = df.dist - np.insert(df.dist.to_numpy(), 0, df.dist[0])[:-1]
        # Adding closer_to_cluster_center column
        if simulated_tempering:
            df['closer_to_cluster_center'] = df.frame.apply(lambda x: True if x in frame_close_to_centroids.values() else False)
        df.to_csv('conf_windows_summary.csv')

        print('Creating the files needed for wham')
        # Creating the files needed for wham.
        if simulated_tempering:
            endung = '_ST'
        else:
            endung = ''
        tpr_files = open(os.path.join(window_path, "tpr_files.dat"), "w")
        pullf_files = open(os.path.join(window_path, "pullf_files.dat"), "w")
        pullx_files = open(os.path.join(window_path, "pullx_files.dat"), "w")
        for window_ID in windows:
            tpr_files.write(f"./{str(window_ID).zfill(5)}/production.tpr\n")
            pullf_files.write(f"./{str(window_ID).zfill(5)}/production_pullf{endung}.xvg\n")
            pullx_files.write(f"./{str(window_ID).zfill(5)}/production_pullx{endung}.xvg\n")
        tpr_files.close()
        pullf_files.close()
        pullx_files.close()

        if simulated_tempering:
            print('Preparing cluster center simulations')
            jobs2launch = {}
            for cluster in frame_close_to_centroids:
                window_ID = frame2window[frame_close_to_centroids[cluster]]
                conf = f"conf{frame_close_to_centroids[cluster]}.gro"
                ID_dir_basename = str(window_ID).zfill(5)
                tools.makedirs(ID_dir_basename)
                ID_dir_abspath = os.path.abspath(ID_dir_basename)
                jobs2launch[ID_dir_abspath] = {'cluster': cluster, 'window_ID':window_ID}
                tools.mv(os.path.join('split_xtc',conf), ID_dir_abspath)
                [tools.cp(new_path_dict[key], ID_dir_abspath, r=True) for key in new_path_dict]

                

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
                    orient_restrain=orient_restrain,
                    orient_restrain_init=orient_restrain_init,
                    orient_restrain_k=orient_restrain_k,
                    orient_restrain_vec=orient_restrain_vec,
                    # General options for the actual pulling
                    pull_nstxout = int(simulation_time['equilibration'] / (dt * number_frame['equilibration'])), # Reduce the output frequency of nstxout. It is not relevant. The frequency was set equal to the nstxout_compressed
                    pull_nstfout = int(simulation_time['equilibration'] / (dt * number_frame['equilibration'])), # Reduce the output frequency of nstfout. It is not relevant. The frequency was set equal to the nstxout_compressed
                    pull_coord1_k=windows[window_ID]['force_constant'],
                    pull_coord1_vec=" ".join([str(xi) for xi in pull_coord_vec]),
                    pull_coord1_init=windows[window_ID]['dist'],
                    pull_coord1_rate= 0.0,)

                # Adding the section for the annealing if required.
                if annealing:
                    # This is a fix method
                    # Increase in 50 K the temperature for the ligands and the protein
                    # during the initial 25 % of the time, then the continued with a cooling (50%)
                    # And finally set constant the temperature to the initial one
                    GroupBool = (np.array(tc_grps) != 'MEMB') & (np.array(tc_grps) != 'SOLV')
                    temp_list = np.arange(temperature, temperature + 52, 2)
                    MDP_equilibration.annealing2(GroupBool, temp_list, heating_frac=0.25, cooling_frac=0.50)
                MDP_equilibration.write(os.path.join(ID_dir_abspath, 'equilibration.mdp'))

                #Create the mdp for get the initial weights
                MDP_annealing_for_ST = mdp.MDP(
                    **mdp_options,
                    define=f"-DPOSRES -DPOSRES_FC_BB=0.0 -DPOSRES_FC_SC=0.0 -DPOSRES_FC_LIPID=0.0 -DDIHRES -DDIHRES_FC=0.0 -DPOSRES_LIG=0.0",
                    time=simulation_time['annealing_for_ST'],
                    xtc_numb_frame=number_frame['annealing_for_ST'])
                MDP_annealing_for_ST.pull(
                    pull_distance,
                    ligands,
                    # Options for the cylinder and orientation restraints
                    flat_bottom_init=flat_bottom_init,
                    flat_bottom_k=flat_bottom_k,
                    orient_restrain=orient_restrain,
                    orient_restrain_init=orient_restrain_init,
                    orient_restrain_k=orient_restrain_k,
                    orient_restrain_vec=orient_restrain_vec,
                    # General options for the actual pulling
                    pull_nstxout = int(simulation_time['annealing_for_ST'] / (dt * number_frame['annealing_for_ST'])), # Reduce the output frequency of nstxout. It is not relevant. The frequency was set equal to the nstxout_compressed
                    pull_nstfout = int(simulation_time['annealing_for_ST'] / (dt * number_frame['annealing_for_ST'])), # Reduce the output frequency of nstfout. It is not relevant. The frequency was set equal to the nstxout_compressed
                    pull_coord1_k=windows[window_ID]['force_constant'],
                    pull_coord1_vec=" ".join([str(xi) for xi in pull_coord_vec]),
                    pull_coord1_init=windows[window_ID]['dist'],
                    pull_coord1_rate= 0.0,)
                # Es importante saber que ecauciones se van a usar. el metodo de pandes con alpha no es del todo correcto dado que lo que se estimo directamente
                # fue el valor de los coeficientes y no de beta F. Con mi ecuacion la formula cambia porque tendriamos que buscar el punto medio no los extremos,
                # ESTO ESTA MAL CONCEPTUALEMNTE PERO LO VOY A USAR PARA PODER SEGUIR PROGRAMANDO LA SIGUIENTEW PARTE
                # The idea is explore a range of temperature from the target tmeperature to 30 K more
                # Por lo pronto lo que puedo hacer es usar ST sin el alpha, que esta bien
                # get the total energies on those points and compute theweights

                if test:
                    annealing_for_ST_temp, annealing_for_ST_time = MDP_annealing_for_ST.annealing3(NumbOfGroups = len(tc_grps), temp = simulated_tempering_temperatures, heating_time=2, constant_temp_time=2)
                else:
                    annealing_for_ST_temp, annealing_for_ST_time = MDP_annealing_for_ST.annealing3(NumbOfGroups = len(tc_grps), temp = simulated_tempering_temperatures, heating_time=100, constant_temp_time=2000)
                MDP_annealing_for_ST.write(os.path.join(ID_dir_abspath, 'annealing_for_ST.mdp'))

                # Create the 1_job.sh
                JOBSH_equilibration_annealing_for_ST = jobsh.JOB(
                    sbatch_keywords={
                        'job-name':f"{window_ID}_getting_weights",
                        'exclude':exclude_nodes,
                        'cpus-per-task': cpu,
                        'partition': partition,
                    },
                    GROMACS_version = GROMACS_version,
                    build_GROMACS_section= f"{str_simulated_tempering_alpha}\n"\
                        f"gmx grompp -f equilibration.mdp -o equilibration.tpr -c {conf} -r {conf} -p {os.path.basename(new_path_dict['topol'])} -n {os.path.basename(new_path_dict['ndx'])}\n"\
                        f"gmx mdrun -nt {cpu} -cpi -stepout 5000 -v -deffnm equilibration -px equilibration_pullx -pf equilibration_pullf {mdrun_opt} {update_gpu} >& equilibration.lis\n\n\n"\
                        f"gmx grompp -f annealing_for_ST.mdp -o annealing_for_ST.tpr -c equilibration.gro -r equilibration.gro -p {os.path.basename(new_path_dict['topol'])} -n {os.path.basename(new_path_dict['ndx'])}\n"\
                        f"gmx mdrun -nt {cpu} -cpi -stepout 5000 -v -deffnm annealing_for_ST -px annealing_for_ST_pullx -pf annealing_for_ST_pullf {mdrun_opt} {update_gpu} >& annealing_for_ST.lis\n"
                )
                JOBSH_equilibration_annealing_for_ST.write(os.path.join(ID_dir_abspath,'1_job.sh'))

            # Launch the 1_job.sh files and wait for end of execution of all center of windows: DONE!
            print('Launching the getting weights jobs (equilibration and annelaing for ST)')
            if frontend:
                tools.job_launch_list([os.path.join(path, '1_job.sh') for path in jobs2launch], shell='bash')# This is just for testing
            else:
                # This will launch and reapeat (jobsh_name) in case that logfile did not finish correctlly due to runing out of tiem on the selected partition,
                process = []
                for job2launch in jobs2launch:
                    process.append(
                        multiprocessing.Process(
                            target=tools.launch_wait_check_repeat,
                            kwargs= {
                                'partition': partition,
                                'jobpaths':[job2launch],
                                'logfile': 'annealing_for_ST.log',
                                'lisfile':'annealing_for_ST.lis',
                                'jobsh_name':'1_job.sh',
                                'logfile_start_datetime': 'equilibration.log',
                            }
                        )
                    )
                [p.start() for p in process]
                [p.join() for p in process]
            print('Finished getting weights jobs')
            # loop over close to centers
            print('Prepearing production simulations for the cluster centers')
            for ID_dir_abspath in jobs2launch:
                energies = weights.get_energy(os.path.join(ID_dir_abspath, 'annealing_for_ST.edr'), annealing_for_ST_time, energy_type = 'Potential', out_fig= os.path.join(ID_dir_abspath,'energy_distribution.svg'))
                # Calculate initial weights for each closest to center
                PandeWeights = weights.PandeWeights(
                    temperatures = sorted(set(annealing_for_ST_temp)),
                    energies = energies,
                    alpha = simulated_tempering_alpha,
                    weight_GROMACS_format = True
                    )
                # Create mdp, assigning initial weights for production and burning phase
                
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
                    orient_restrain=orient_restrain,
                    orient_restrain_init=orient_restrain_init,
                    orient_restrain_k=orient_restrain_k,
                    orient_restrain_vec=orient_restrain_vec,
                    # General options for the actual pulling
                    pull_coord1_k=windows[window_ID]['force_constant'],
                    pull_coord1_vec=" ".join([str(xi) for xi in pull_coord_vec]),
                    pull_coord1_init=windows[window_ID]['dist'],
                    pull_coord1_rate= 0.0,)
                MDP_production.simulated_tempering(
                    temp_list = simulated_tempering_temperatures, # The same used for the anealing for the prediction of the inital weights
                    init_lambda_weights = PandeWeights,
                )
                MDP_production.write(os.path.join(ID_dir_abspath,'production.mdp'))
                # Create the 2_job.sh, GROMACS only support -update gpu for the integtration method 'md' and becasue we are doing ST ('md-vv'), we are not able to use it. This is a lack of efficience.
                JOBSH_production = jobsh.JOB(
                    sbatch_keywords={
                        'job-name':f"{jobs2launch[ID_dir_abspath]['window_ID']}_production",
                        'exclude':exclude_nodes,
                        'cpus-per-task': cpu,
                        'partition': partition,
                    },
                    GROMACS_version = GROMACS_version,
                    build_GROMACS_section=f"{str_simulated_tempering_alpha}\n"\
                        f"gmx grompp -f production.mdp -o production.tpr -c equilibration.gro -r equilibration.gro -p {os.path.basename(new_path_dict['topol'])} -n {os.path.basename(new_path_dict['ndx'])}\n"\
                        f"gmx mdrun -nt {cpu} -cpi -stepout 5000 -v -deffnm production -px production_pullx -pf production_pullf {mdrun_opt} >& production.lis\n"
                )
                JOBSH_production.write(os.path.join(ID_dir_abspath,'2_job.sh'))
            print('Finished preparetion of production simulations for the cluster centers')
            # For this part multithreading will be used in order to spend less time wating
            # print(f'Creating the threads. Number of threads before opening: {threading.activeCount()}')
            if frontend:
                # This is just for testing
                center_process = [
                    multiprocessing.Process(
                        target=tools.job_launch_list,
                        kwargs = {
                            'job_path_list':[os.path.join(path, '2_job.sh') for path in jobs2launch],
                            'shell':'bash'
                        }
                    )
                ]

            else:
                # This will launch and reapeat (jobsh_name) in case that logfile did not finish correctlly due to runing out of time on the selected partition,
                center_process = []
                for job2launch in jobs2launch:
                    center_process.append(
                        multiprocessing.Process(
                            target=tools.launch_wait_check_repeat,
                            kwargs= {
                                'partition': partition,
                                'jobpaths':[job2launch],
                                'logfile': 'production.log',
                                'lisfile':'production.lis',
                                'jobsh_name':'2_job.sh',
                            }
                        )
                    )
            # print(f'The main thread was set. This is the launch of the production simulations on the cluster centers which corresponds to launch_wait_check_repeat() in case of frontend True or job_launch_list() if False . Number of threads right now (after set but not start): {threading.activeCount()}')
            member_process = []
            # Iterate over the cluster center (jobs2launch)
            for ID_dir_abspath in jobs2launch:
                # Create the Procese for every cluster
                # Iterate over the member of the corresponded cluster center
                for frame in user_membership[jobs2launch[ID_dir_abspath]['cluster']]:
                    setup_production_kwargs = {
                        'frame':frame,
                        'frame2window':frame2window,
                        'log_center_cluster_abspath':os.path.join(ID_dir_abspath, 'production.log'),
                        'simulated_tempering_weights_equilibration_time':simulated_tempering_weights_equilibration_time,
                        'new_path_dict':new_path_dict,
                        'mdp_options':mdp_options,
                        'simulation_time':simulation_time,
                        'number_frame':number_frame,
                        'pull_distance':pull_distance,
                        'ligands':ligands,
                        'flat_bottom_init':flat_bottom_init,
                        'flat_bottom_k':flat_bottom_k,
                        'orient_restrain': orient_restrain,
                        'orient_restrain_init':orient_restrain_init,
                        'orient_restrain_k':orient_restrain_k,
                        'orient_restrain_vec':orient_restrain_vec,
                        'dt':dt,
                        'windows':windows,
                        'pull_coord_vec':pull_coord_vec,
                        'annealing':annealing,
                        'tc_grps':tc_grps,
                        'temperature':temperature,
                        'exclude_nodes':exclude_nodes,
                        'cpu':cpu,
                        'partition':partition,
                        'GROMACS_version':GROMACS_version,
                        'str_simulated_tempering_alpha':str_simulated_tempering_alpha,
                        'simulated_tempering_temperatures':simulated_tempering_temperatures,
                        'frontend':frontend,
                        'mdrun_opt':mdrun_opt,
                    }
                    # Creating as many process as windows are in the simulation, the identificator is the window ID
                    member_process.append(multiprocessing.Process(target = setup_production, kwargs=setup_production_kwargs))
            # Set what kind of launch we need if frontend or not is used. In the first case is just a normal synchronic job. The secondone is asynchronic jobs.
            # print(f'Starting the threads and doing the corresponded join sequence. Number of threads: {threading.activeCount()}')
            if frontend:
                for center_p in center_process:
                    center_p.start()
                    center_p.join()
                for member_p in member_process:
                    member_p.start()
                    member_p.join()
            else:
                # Start the threads
                [center_p.start() for center_p in center_process]
                [member_p.start() for member_p in member_process]
                # print(f'Number of theads before the join (after the start): {threading.activeCount()}')
                [center_p.join() for center_p in center_process]
                [member_p.join() for member_p in member_process]
                # print(f'Number of theads after the join: {threading.activeCount()}')

            # print(f'Number of theads when everything was finished): {threading.activeCount()}')
            print('Everything finished, now the analysis.')

            for window_ID in windows:
                print(f'Analysing window: {window_ID}')
                ID_dir_abspath = os.path.abspath(str(window_ID).zfill(5))
                production_xvg_data = xvg.XVG(os.path.join(ID_dir_abspath, 'production.xvg')).data
                lowest_temperature_indexes = production_xvg_data[:,1] == 0 # Those points that belongs to state 0
                # From production_xvg I take the info for the construction of the Markov Chain matrix.
                # I need to get the frequency of output (nstdhdl) and the frequency of attempted moves changing the temperature (nstexpanded).
                # For it, we take this parameters from the last MDP_production (all are the same)
                nstdhdl = MDP_production.keywords['Simulated Tempering']['nstdhdl']
                nstexpanded = MDP_production.keywords['Simulated Tempering']['nstexpanded']
                if nstexpanded > nstdhdl:
                    if nstdhdl % nstexpanded:
                        markov.MARKOV(production_xvg_data[::int(nstexpanded/nstdhdl),1].astype(int)).plot(os.path.join(ID_dir_abspath, 'transition_matrix.svg'))
                    else:
                        print('The transition matrix were not export because nstexpanded is not multiple of nstdhdl.')
                elif nstdhdl == nstexpanded:
                    markov.MARKOV(production_xvg_data[:,1].astype(int)).plot(os.path.join(ID_dir_abspath, 'transition_matrix.svg'))
                else:
                    print('The transition matrix were not export because nstexpanded < nstdhdl.')

                # Plot the weights convergence
                _, _ = weights.get_weights_from_log(os.path.join(ID_dir_abspath,'production.log'), plot = True)

                # Taking only the data for the lowest temperature
                pullf_xvg = xvg.XVG(os.path.join(ID_dir_abspath,'production_pullf.xvg'))
                pullf_xvg.data = pullf_xvg.data[lowest_temperature_indexes]
                pullf_xvg.write(os.path.join(ID_dir_abspath,'production_pullf_ST.xvg'))

                pullx_xvg = xvg.XVG(os.path.join(ID_dir_abspath,'production_pullx.xvg'))
                pullx_xvg.data = pullx_xvg.data[lowest_temperature_indexes]
                pullx_xvg.write(os.path.join(ID_dir_abspath,'production_pullx_ST.xvg'))
            print('The analysis finished')
        else:
            jobs2launch = []
            for window_ID in windows:
                frame = windows[window_ID]['frame']
                conf = f"conf{frame}.gro"
                ID_dir_basename = str(window_ID).zfill(5)
                tools.makedirs(ID_dir_basename)
                ID_dir_abspath = os.path.abspath(ID_dir_basename)
                jobs2launch.append(ID_dir_abspath)
                tools.mv(os.path.join('split_xtc',conf), ID_dir_abspath)
                [tools.cp(new_path_dict[key], ID_dir_abspath, r=True) for key in new_path_dict]

                

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
                    orient_restrain=orient_restrain,
                    orient_restrain_init=orient_restrain_init,
                    orient_restrain_k=orient_restrain_k,
                    orient_restrain_vec=orient_restrain_vec,
                    # General options for the actual pulling
                    pull_nstxout = int(simulation_time['equilibration'] / (dt * number_frame['equilibration'])), # Reduce the output frequency of nstxout. It is not relevant. The frequency was set equal to the nstxout_compressed
                    pull_nstfout = int(simulation_time['equilibration'] / (dt * number_frame['equilibration'])), # Reduce the output frequency of nstfout. It is not relevant. The frequency was set equal to the nstxout_compressed
                    pull_coord1_k=windows[window_ID]['force_constant'],
                    pull_coord1_vec=" ".join([str(xi) for xi in pull_coord_vec]),
                    pull_coord1_init=windows[window_ID]['dist'],
                    pull_coord1_rate= 0.0,)

                # Adding the section for the annealing if required.
                if annealing:
                    # This is a fix method
                    # Increase in 50 K the temperature for the ligands and the protein
                    # during the initial 25 % of the time, then the continued with a cooling (50%)
                    # And finally set constant the temperature to the initial one
                    GroupBool = (np.array(tc_grps) != 'MEMB') & (np.array(tc_grps) != 'SOLV')
                    temp_list = np.arange(temperature, temperature + 52, 2)
                    MDP_equilibration.annealing2(GroupBool, temp_list, heating_frac=0.25, cooling_frac=0.50)

                MDP_equilibration.write(os.path.join(ID_dir_abspath, 'equilibration.mdp'))

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
                    orient_restrain=orient_restrain,
                    orient_restrain_init=orient_restrain_init,
                    orient_restrain_k=orient_restrain_k,
                    orient_restrain_vec=orient_restrain_vec,
                    # General options for the actual pulling
                    pull_coord1_k=windows[window_ID]['force_constant'],
                    pull_coord1_vec=" ".join([str(xi) for xi in pull_coord_vec]),
                    pull_coord1_init=windows[window_ID]['dist'],
                    pull_coord1_rate= 0.0,)

                MDP_production.write(os.path.join(ID_dir_abspath, 'production.mdp'))

                # I will need to modify templates.GROMACS_WINDOW_JOB and adding the posibility of annealing caclulation
                # Also, I will not be able to lanch as I do all the jobs at once becasue some jobs are need to be launch first
                # And then the rest one
                JOBSH_window = jobsh.JOB(
                    sbatch_keywords={
                        'job-name':f"{window_ID}_production",
                        'exclude':exclude_nodes,
                        'cpus-per-task': cpu,
                        'partition': partition,
                    },
                    GROMACS_version = GROMACS_version,
                    build_GROMACS_section=f"{str_simulated_tempering_alpha}\n"\
                        f"gmx grompp -f equilibration.mdp -o equilibration.tpr -c {conf} -r {conf} -p {os.path.basename(new_path_dict['topol'])} -n {os.path.basename(new_path_dict['ndx'])}\n"\
                        f"gmx mdrun -nt {cpu} -cpi -stepout 5000 -v -deffnm equilibration -px equilibration_pullx -pf equilibration_pullf {mdrun_opt} {update_gpu} >& equilibration.lis\n\n\n"\
                        f"gmx grompp -f production.mdp -o production.tpr -c equilibration.gro -r equilibration.gro -p {os.path.basename(new_path_dict['topol'])} -n {os.path.basename(new_path_dict['ndx'])}\n"\
                        f"gmx mdrun -nt {cpu} -cpi -stepout 5000 -v -deffnm production -px production_pullx -pf production_pullf {mdrun_opt} {update_gpu} >& production.lis\n"
                )
                JOBSH_window.write(os.path.join(ID_dir_abspath, 'job.sh'))

            if frontend:
                tools.job_launch(shell='bash')  # This is just for testing
            else:
                process = []
                for job2launch in jobs2launch:
                    process.append(
                        multiprocessing.Process(
                            target=tools.launch_wait_check_repeat,
                            kwargs= {
                                'partition': partition,
                                'jobpaths':[job2launch],
                                'logfile': 'production.log',
                                'lisfile':'production.lis',
                                'jobsh_name':'job.sh',
                                'logfile_start_datetime': 'equilibration.log',
                            }
                        )
                    )
                [p.start() for p in process]
                [p.join() for p in process]

        # Printing the coordanates.
        number_windows = len(windows)
        # THis case I use 3 because I use 2 pull restrain peer coord. orientation and loss of channel entrance
        if orient_restrain:
            number_coords = 3*len(ligands)
        else:
            number_coords = 2*len(ligands)
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

        # wham_individual_coords = [f"gmx wham -ac -temp {temperature} -zprof0 {pull_distance} -bins 300 -unit kJ -nBootstrap 200 -bs-method hist -is {s} -ix pullx_files.dat -it tpr_files.dat -o {s.split('.')[0]} -hist hist_{s.split('.')[0]} -oiact iact_{s.split('.')[0]}.xvg -bsres bsResult_{s.split('.')[0]}.xvg -bsprof bsProfs_{s.split('.')[0]}.xvg" for s in select_file_names]
        # This is jut to see the histogrma
        tools.run(f"gmx wham -temp {temperature} -bins 300 -unit kJ -is coord0_selected.dat -ix pullx_files.dat -it tpr_files.dat -hist hist_coord0_selected -histonly -b 0")

    # !!!! Activate an automatic function that zip/tar all the files, except the xtc, but I need to conserv the directory tree.
    os.chdir(cwd)
    end_datetime = datetime.datetime.now()
    print(
        f"The function automation.main() finished. The script start at {start_datetime} and ended at {end_datetime}; the execution time was: {round((end_datetime - start_datetime).total_seconds() / 86400, 3)} days!")


if __name__ == '__main__':
    """This is just the main function of the script executing top_parser, atom_block_modifier
    and writer functions, see __doc__ for details"""
    pass

    # parser = argparse.ArgumentParser(description=__doc__,
    #                                  formatter_class=argparse.RawTextHelpFormatter)

    # parser.add_argument('--input_path_dict',
    #                     help="String-like-dict. This dictionary is refereed to the one that you get from: dict_gen()['receptor_name']['ligand_name']",
    #                     dest='input_path_dict',
    #                     type=str)
    # parser.add_argument('--output_path',
    #                     help="Output path",
    #                     dest='output_path',
    #                     type=str)
    # parser.add_argument('--lig_vec_group1',
    #                     help="The index atoms of Ligand to consider. This are could be taken from a pymol session the vector will e constructed from 1-->2",
    #                     nargs="+",
    #                     dest='lig_vec_group1',
    #                     type=int)
    # parser.add_argument('--lig_vec_group2',
    #                     help="The index atoms of Ligand to consider. This are could be taken from a pymol session the vector will e constructed from 1-->2",
    #                     nargs="+",
    #                     dest='lig_vec_group2',
    #                     type=int)
    # parser.add_argument('--dt',
    #                     help="Integration time in ps",
    #                     dest='dt',
    #                     default=0.004,
    #                     type=float)
    # parser.add_argument('--simulation_time',
    #                     help="There are three different simulation: pulling, window equilibration, window production."
    #                     "You must provied the three corresponded times in this order (in ns). Defaults are:"
    #                     "pulling = 100000 ps"
    #                     "equilibration = 2000 ps"
    #                     "production = 20000 ps",
    #                     dest='simulation_time',
    #                     nargs=3,
    #                     default=[100000, 2000, 20000],
    #                     type=float)
    # parser.add_argument('--number_frame',
    #                     help="There are three different simulation: pulling, window equilibration, window production."
    #                     "You must provied the three corresponded number of frames in this order. Defaults are:"
    #                     "pulling = 2000"
    #                     "equilibration = 2"
    #                     "production = 2000",
    #                     dest='number_frame',
    #                     nargs=3,
    #                     default=[2000, 2, 2000],
    #                     type=int)
    # parser.add_argument('--temperature',
    #                     help="Temperature in K",
    #                     dest='temperature',
    #                     default=303.15,
    #                     type=float)
    # parser.add_argument('--pull_coord_vec',
    #                     help="This is the vector along the pulling will be done. Default is [0.0, 0.0, -1.0]",
    #                     dest='pull_coord_vec',
    #                     nargs=3,
    #                     default=[0.0, 0.0, -1.0],
    #                     type=float)
    # parser.add_argument('--flat_bottom_init',
    #                     help="The radium for the flat bottom restrain. Default is 0.5 nm",
    #                     dest='flat_bottom_init',
    #                     default=0.5,
    #                     type=float)
    # parser.add_argument('--flat_bottom_k',
    #                     help="The force constant for the flat bottom restrain. Default is 400 kJ/mol*nm^2",
    #                     dest='flat_bottom_k',
    #                     default=400,
    #                     type=float)
    # parser.add_argument('--orient_restrain_init',
    #                     help="The angle for the orientation restrain. Default is 0.349 rad, 20 degrees 20/180*3.141",
    #                     dest='orient_restrain_init',
    #                     default=0.349,
    #                     type=float)
    # parser.add_argument('--orient_restrain_k',
    #                     help="The force constant for the orientation restrain. Default is 1308.36 kJ/mol/rad^2, leads to sigma = 2.5 degrees",
    #                     dest='orient_restrain_k',
    #                     default=1308.36,
    #                     type=float)
    # parser.add_argument('--orient_restrain_vec',
    #                     help="This means that the vector will be calculate in suing the initial orientation fo the ligands. if a vector is provied, for example (0,0,-1), then this will be used. You need to provied the vector in the format '0 0 -1'. 'calculate' is default",
    #                     dest='orient_restrain_vec',
    #                     default="calculate",
    #                     type=str)
    # parser.add_argument('--pull_distance',
    #                     help="''",
    #                     dest='pull_distance',
    #                     default=4.5,
    #                     type=float)
    # parser.add_argument('--pull_middle_points',
    #                     help="''",
    #                     nargs="+",
    #                     dest='pull_middle_points',
    #                     default=[],
    #                     type=float)
    # parser.add_argument('--pull_force_constants',
    #                     help="''",
    #                     nargs="+",
    #                     dest='pull_force_constants',
    #                     default=[10000],
    #                     type=float)
    # parser.add_argument('--window_widths',
    #                     help="''",
    #                     nargs="+",
    #                     dest='window_widths',
    #                     default=[0.1],
    #                     type=float)
    # parser.add_argument('--chains',
    #                     help="''",
    #                     nargs="+",
    #                     dest='chains',
    #                     default=['A', 'B', 'C', 'D', 'E'],
    #                     type=str)
    # parser.add_argument('--GROMACS_version',
    #                     help="",
    #                     dest='GROMACS_version',
    #                     default="2021.1",
    #                     type=str)
    # parser.add_argument('--cpu',
    #                     help="",
    #                     dest='cpu',
    #                     default=12,
    #                     type=int)
    # parser.add_argument('--COM_dist_cpu',
    #                     help="",
    #                     dest='COM_dist_cpu',
    #                     default=12,
    #                     type=int)
    # parser.add_argument('--frontend',
    #                     help='Type of machine smaug, gwdg or frontend',
    #                     dest='frontend',
    #                     default=False,
    #                     type=bool)
    # parser.add_argument('--exclude_nodes',
    #                     help="",
    #                     dest='exclude_nodes',
    #                     default="",
    #                     type=str)
    # parser.add_argument('--pull_code',
    #                     help='If you want to execute the pulling code add the flag . Default is False',
    #                     nargs="?",
    #                     dest='pull_code',
    #                     const=True,
    #                     default=False,
    #                     type=bool)
    # parser.add_argument('--window_code',
    #                     help='If you want to execute the window code. Default is False',
    #                     nargs="?",
    #                     dest='window_code',
    #                     const=True,
    #                     default=False,
    #                     type=bool)
    # args = parser.parse_args()

    # args.input_path_dict = eval(args.input_path_dict)
    # args.simulation_time = {
    #     'pulling': args.simulation_time[0], 'equilibration': args.simulation_time[1], 'production': args.simulation_time[2]}
    # args.number_frame = {
    #     'pulling': args.number_frame[0], 'equilibration': args.number_frame[1], 'production': args.number_frame[2]}
    # args.orient_restrain_vec = args.orient_restrain_vec.strip().lower()
    # if args.orient_restrain_vec != "calculate":
    #     args.orient_restrain_vec = [
    #         float(xi) for xi in args.orient_restrain_vec.split()]

    # main(args.input_path_dict,
    #      args.output_path,
    #      args.lig_vec_group1,
    #      args.lig_vec_group2,
    #      dt=args.dt,
    #      simulation_time=args.simulation_time,
    #      number_frame=args.number_frame,
    #      temperature=args.temperature,
    #      pull_coord_vec = args.pull_coord_vec,
    #      flat_bottom_init=args.flat_bottom_init,
    #      flat_bottom_k=args.flat_bottom_k,
    #      orient_restrain_init=args.orient_restrain_init,
    #      orient_restrain_k=args.orient_restrain_k,
    #      orient_restrain_vec=args.orient_restrain_vec,
    #      pull_distance=args.pull_distance,
    #      pull_middle_points=args.pull_middle_points,
    #      pull_force_constants=args.pull_force_constants,
    #      window_widths=args.window_widths,
    #      chains=args.chains,
    #      GROMACS_version=args.GROMACS_version,
    #      cpu=args.cpu,
    #      COM_dist_cpu=args.COM_dist_cpu,
    #      frontend=args.frontend,
    #      exclude_nodes=args.exclude_nodes,
    #      pull_code=args.pull_code,
    #      window_code=args.window_code)
