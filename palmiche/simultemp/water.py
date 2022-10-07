from palmiche.tools import tools, mdp, xvg, markov
from palmiche.simultemp import weights
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import socket


def main(wd = './',
        weights_method = 'PandeWeights', #'AleWeights', 'PandeWeights'
        EnergyType = 'Potential', # 'Total-Energy'
        heat_fraction = 0.25, # DEpendding on what weights_method is chosen heat_fraction or heating_time and constant_temp_time will be used
        heating_time = 100,
        constant_temp_time = 2000,
        number_of_tc_grps = 1,
        temp_list = range(300, 332, 2),
        maxsol = 500,
        box = (2.5,2.5,2.5),
        water_ff = 'tip3p',
        ions_ff = 'oplsaa',
        simulation_time = {'annealing':60000,'tempering': 60000}, # ps the annealing time only will be used for the method AleWeights, and is a similar time to the one use for the calculation of PandeWeights
        xtc_numb_frame = 1000,
        weight_GROMACS_format = True,
        cpu = 12,
        # This arguments will give the option to only collect those point in an specific 
        # interval of time. E. g = (500, 670) will collect those point between 500 and 670 ps,
        # None means that all point will be collected. This feature could be useful in case that
        # only some region of the simulation want to be analyzed (the coling region, for example).
        # the interval is closed, both extrames are taken
        collecting_interval = None,
        SIMTEMP_OCCUP_STATE0 = 0,
        **mdp_options,):

    if SIMTEMP_OCCUP_STATE0:
       str_SIMTEMP_OCCUP_STATE0 = f"""
       source /data/shared/opt/gromacs/2021.6_simTemp_setOccupState0/bin/GMXRC.bash
       PATH=$PATH:"/data/shared/opt/gromacs/2021.6_simTemp_setOccupState0/bin"
       export SIMTEMP_OCCUP_STATE0={SIMTEMP_OCCUP_STATE0}
       """
    else:
        str_SIMTEMP_OCCUP_STATE0 = '#' 

    tools.makedirs(wd)
    hostname = socket.gethostname()
    if hostname == 'smaug':
        mdrun_opt = '-pin on -pinstride 1 -pinoffset 0 -gpu_id 0 -nice 19'
    else:
        mdrun_opt = '-pin on -pinstride 1 -pinoffset 0 -nice 19'

    minimization_MDP =mdp.MDP(type = 'minimization',**mdp_options)
    minimization_MDP.write(os.path.join(wd,'minimization.mdp'))

    annealing_MDP = mdp.MDP(
        type = 'production',
        time = simulation_time['annealing'],
        xtc_numb_frame = xtc_numb_frame,
        **mdp_options)
    # Choose what kind of annealing strategy depending the weight calculation method
    if weights_method == 'PandeWeights':
        annealing_temp, annealing_time = annealing_MDP.annealing3(number_of_tc_grps,temp_list, heating_time = heating_time, constant_temp_time = constant_temp_time)
    elif weights_method == 'AleWeights':
        annealing_MDP.annealing1(number_of_tc_grps,temp_list,heat_fraction=heat_fraction,  nstenergy_same_as_nstcalcenergy = True)
    annealing_MDP.write(os.path.join(wd,'annealing.mdp'))

    cwd = os.getcwd()
    os.chdir(wd)

    tools.run(f"""
        {str_SIMTEMP_OCCUP_STATE0}
        export GMX_MAXBACKUP=-1
        gmx solvate -box {' '.join([str(b) for b in box])} -maxsol {maxsol} -o water.gro
        gmx pdb2gmx -f water.gro -o water.gro -water {water_ff} -ff {ions_ff} -p topol.top
        echo q | gmx make_ndx -f water.gro -o index.ndx

        # This is for testing, it should be at the job.sh, and them add the wait function of tools module
        # Or even better launch the script itself. And them everything is done nicely. No need to wait an all of that
        gmx grompp -f minimization.mdp  -c water.gro -p topol.top -n index.ndx -o minimization.tpr
        gmx mdrun -nt {cpu} -stepout 5000 -v -deffnm minimization -cpi -maxh 48 {mdrun_opt} >& minimization.lis
        gmx grompp -f annealing.mdp  -c minimization.gro -p topol.top -n index.ndx -o annealing.tpr
        gmx mdrun -nt {cpu} -stepout 5000 -v -deffnm annealing -cpi -maxh 48 {mdrun_opt} >& annealing.lis
        """)
    
    if weights_method == 'PandeWeights':
        energy = weights.get_energy('annealing.edr', annealing_time, energy_type = EnergyType)

        W = weights.PandeWeights(
                temperatures = sorted(set(annealing_temp)),
                energies = energy,
                alpha = SIMTEMP_OCCUP_STATE0, 
                weight_GROMACS_format = weight_GROMACS_format
                )

    elif weights_method == 'AleWeights':

        """
        Getting the function E(T), the way that will be done is averaging the internal energy 
        The same value of temperature could have several values of total energy
        The first an easiest approximation is take all the simulation into account
        The second one is just take into account the cooling region
        """

        tools.run(f"""
        {str_SIMTEMP_OCCUP_STATE0}
        export GMX_MAXBACKUP=-1
        echo Temperature | gmx energy -f annealing.edr -s annealing.tpr -o temperature.xvg
        echo {EnergyType} | gmx energy -f annealing.edr -s annealing.tpr -o energy.xvg
        """)
        temperature = xvg.XVG('temperature.xvg').data
        energy = xvg.XVG('energy.xvg').data

        # In case that collecting_interval is provided
        if collecting_interval:
            collecting_interval = np.array(collecting_interval)*1000
            # Could also been selected energy because the first column (time) is the same
            mask = (temperature[:,0] >= collecting_interval.min()) & (temperature[:,0] <= collecting_interval.max())
            temperature = temperature[mask]
            energy = energy[mask]


        # df = pd.DataFrame(np.vstack([temperature[:,1], energy[:,1]]).T, columns = ['temperature', 'energy'])
        # mean = df.groupby(['temperature']).mean()
        # mean.sort_index(inplace=True)
        # # # To get a column named temperature
        # # mean['temperature'] = mean.index
        # # mean.reset_index(drop=True, inplace=True)
    
        # Only take into account the study interval of temperature
        # The "-2" rectification is to take into account a little more at the edges, but it will not affect the calculation of the weights
        # It will only add extra points. The problem is that with the filter let say that we add all the number greater than or equal to 300 K
        # But then in the integration 300 is not in the DataFrame, is just for precaution.
        
        mask = (temperature[:,1] >= min(temp_list) - 2) & (temperature[:,1] <= max(temp_list) + 2)
        temperature = temperature[mask]
        energy = energy[mask]

        fig, ax = plt.subplots(figsize = (25,9))
        ax.scatter(temperature[:,1], energy[:,1])
        ax.set(
            xlim = [min(temp_list), max(temp_list)],
            xlabel = 'T [K]',
            ylabel = f'{EnergyType} [kJ/mol]',
            title = f'Dependency of the {EnergyType} with the Temperature'
        )
        fig.savefig(f'{EnergyType}_vs_T.png')
        
        # Calculating the weights
        W = weights.AleWeights(
            temperature[:,1],
            energy[:,1],
            temp_list,
            alpha = SIMTEMP_OCCUP_STATE0,
            weight_GROMACS_format = weight_GROMACS_format)

    # Simulated Tempering
    tempering_MDP = mdp.MDP(
        type = 'production',
        time = simulation_time['tempering'],
        xtc_numb_frame = xtc_numb_frame,
        **mdp_options)
    tempering_MDP.simulated_tempering(temp_list, W)
    tempering_MDP.write(os.path.join('tempering.mdp'))

    tools.run(f"""
        {str_SIMTEMP_OCCUP_STATE0}
        gmx grompp -f tempering.mdp  -c minimization.gro -p topol.top -n index.ndx -o tempering.tpr
        gmx mdrun -nt {cpu} -stepout 5000 -v -deffnm tempering -cpi -maxh 48 {mdrun_opt} >& tempering.lis
        """)
    weights.get_weights_from_log('tempering.log', plot = True)
    
    tempering_xvg_data = xvg.XVG('tempering.xvg').data
    markov.MARKOV(tempering_xvg_data[:,1].astype(int)).plot('transition_matrix.svg')


    os.chdir(cwd)

if __name__ == '__main__':
    #ant the value to 20000 to get more resolution on the energy on this step
    # I need to improve this part because n the simulated tempering I dont need such resolution on the energy
    tc_grps = ['Water']
    ref_t = 300
    # I put const
    mdp_options = {             
    "tc_grps":' '.join(tc_grps),
    "tau_t":' '.join(len(tc_grps)*['1.0']),
    "ref_t":' '.join(len(tc_grps)*[str(ref_t)]),
    "dt": 0.002,
    "pcoupltype": 'isotropic',
    "compressibility": 4.5e-5,
    "ref_p": 1.0,
    #'rlist': 1,
    #'rcoulomb': 1,
    }
    args = [
            {
                    'wd' : 'PandeWeights',
                    'weights_method' : 'PandeWeights',
                    'EnergyType':'Potential', # 'Total-Energy'
                    'heat_fraction' : 0.25,
                    'heating_time' : 2,
                    'constant_temp_time' : 5,
                    'number_of_tc_grps' : 1,
                    'temp_list' : range(300, 332, 2),
                    'maxsol' : 500,
                    'box' : (2.5,2.5,2.5),
                    'water_ff' : 'tip3p',
                    'ions_ff' : 'oplsaa',
                    'simulation_time' : {'annealing':300,'tempering': 50}, #ps
                    'xtc_numb_frame' : 1000,
                    'weight_GROMACS_format' : True,
                    'cpu' : 12,
                    'collecting_interval' : None,
                    'SIMTEMP_OCCUP_STATE0': 0,
            },
            {
                    'wd' : 'AleWeights',
                    'weights_method' : 'AleWeights',
                    'EnergyType':'Potential', # 'Total-Energy'
                    'heat_fraction' : 0.25,
                    'heating_time' : 2,
                    'constant_temp_time' : 5,
                    'number_of_tc_grps' : 1,
                    'temp_list' : range(300, 332, 2),
                    'maxsol' : 500,
                    'box' : (2.5,2.5,2.5),
                    'water_ff' : 'tip3p',
                    'ions_ff' : 'oplsaa',
                    'simulation_time' : {'annealing':300,'tempering': 50}, #ps
                    'xtc_numb_frame' : 1000,
                    'weight_GROMACS_format' : True,
                    'cpu' : 12,
                    'collecting_interval' : None,    
                    'SIMTEMP_OCCUP_STATE0':0, 
            }
    ]
    for arg in args:
        arg.update(mdp_options)
        main(**arg)