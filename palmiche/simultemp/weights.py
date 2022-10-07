#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from palmiche.utils import tools, xvg
import tempfile
import os
import tqdm
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import integrate
import warnings

def AleWeights(temperatures, internal_energies, target_temperatures, integrator = 'trapezoid', alpha = None, weight_GROMACS_format = True):
    """
    Here explain the definition of the gromacs weight
    https://manual.gromacs.org/documentation/2019/reference-manual/algorithms/expanded-ensemble.html
    Tengo que mejorar lo del input del data frame es muy complicado dar ese input, que sea una arraytemperatura, el otro el valor de la energia, y em temp_list
    para las temperaturasw que quiero calcular los weights, de esta forma es mucho mas simple, tengo que ver el signo de la matriz porque sino tengo que cambiar el signo par lo del alpha
    MEJORAR!!!!
    Calculation of the weights for simulated tempering. The weights respond to the hamiltonian:
        H(p,q,T) = V(q) + K(p) - weight/beta; where beta = 1/KbT
    In GROMACS the hamiltonian is deffined:
        H(p,q,T) = V(q) + K(p) - weight_GROMACS
    Therefore:
        weight_GROMACS = weight/beta
    To convert to GROMACS weights you just need to provided weight_GROMACS_format = True (default option)

    Args:
        DataFrame (pd.DataFrame): A Data Frame with temperatures as index and one column with the internal energy with name total_energy.
        temp_list (list): The list of temperature to calculate the weights. At least two. If some temperature is not in the temperature range provided with DataFrame a ValueError will prompt.
        integrator (str, optional): Integration method. Are currently supported:
            trapezoid --> scipy.integrate.trapezoid
            trapz --> numpy.trapz (Default)
            simps --> scipy.integrate.simpson 
        weight_GROMACS_format (bool, optional): If True each weight will be divided by its corresponded beta value of temperature. Defaults to True.

    Returns:
        weights (numpy array): weights for each temperature in the order of temp_list. weights[0] will be the weight of temp_list[0] and so on.
    """

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
    target_temperatures = np.array(target_temperatures)
    # Converting to pandas.DataFrame and sorting the values respect to temperature
    DataFrame = pd.DataFrame(
        {
           'temperatures': temperatures,
           'internal_energies': internal_energies,
        }
    )
    DataFrame.sort_values(by='temperatures', inplace=True)
    DataFrame.reset_index(drop=True, inplace=True)
    # Converting target_temperatures to numpy array
    # Creating the integrand U(T)/T^2
    DataFrame['integrand'] = DataFrame.internal_energies / (DataFrame.temperatures**2)
    T_MAX = DataFrame.temperatures.max()
    T_MIN = DataFrame.temperatures.min()


    # Deffining the coefficient matrix
    M = 1*np.eye(len(target_temperatures), k = 1, dtype=int)
    np.fill_diagonal(M, -1)
    M[-1] = np.ones(len(target_temperatures), dtype = int)

    # Calculating solution vector, based on integral from Ta to Tb of 1/KB * U(T)/T^2
    # Calculating the first component depending if alpha was selected or not
    print('Calculating integrals for solution vector...')
    T0, T1 = target_temperatures[0], target_temperatures[1]
    if not T_MIN <= T0 <= T_MAX:
        warnings.warn(f"Temperature {T0} provided with temp_list is not in the interval of DataFrame; which is [{T_MIN}, {T_MAX}]")
    if not T_MIN <= T0 <= T_MAX:
        warnings.warn(f"Temperature {T1} provided with temp_list is not in the interval of DataFrame; which is [{T_MIN}, {T_MAX}]")
    # Selecting the interval of integration
    data = DataFrame[(DataFrame.temperatures >= T0) & (DataFrame.temperatures <= T1)]
    if len(data) < 4: # At least 4 points
        raise RuntimeError(f"There are not enough data points for the range of temperatures {T0} to {T1}")

    if alpha:
        gamma = alpha/(1-alpha)
        b = [-integrator(data.integrand, data.temperatures) / tools.CTE.Kb - np.log(gamma*(len(temperatures) - 1))]
    else:
        b = [-integrator(data.integrand, data.temperatures) / tools.CTE.Kb]

    for i in tqdm.tqdm(range(1, len(target_temperatures)), total =  len(target_temperatures) - 2):
        if i < len(target_temperatures) - 1:
            # Getting integration interval extrem values and checking that belongs to the DataFrame provided.
            T0, T1 = target_temperatures[i], target_temperatures[i+1]
            if not T_MIN <= T0 <= T_MAX:
                warnings.warn(f"Temperature {T0} provided with temp_list is not in the interval of DataFrame; which is [{T_MIN}, {T_MAX}]")
            if not T_MIN <= T0 <= T_MAX:
                warnings.warn(f"Temperature {T1} provided with temp_list is not in the interval of DataFrame; which is [{T_MIN}, {T_MAX}]")

            # Selecting the interval of integration
            data = DataFrame[(DataFrame.temperatures >= T0) & (DataFrame.temperatures <= T1)]
            if len(data) < 4: # At least 4 points
                raise RuntimeError(f"There are not enough data points for the range of temperatures {T0} to {T1}")

            b.append(-integrator(data.integrand, data.temperatures) / tools.CTE.Kb)
        else:
            # We reach the last temperature
            b.append(0)
            break
    # Calculate the weights:
    weights = np.linalg.solve(M, b)
     
    # Convert if needed to GROMACS weights
    if weight_GROMACS_format:
        minimum = weights.min()
        if minimum > 0:
            weights -= minimum
        else:
            weights += np.abs(minimum)
    
    return weights

def PandeWeights(temperatures, energies, alpha = None, weight_GROMACS_format = True):
    # Here explain the definition of the gromacs weight
    # https://manual.gromacs.org/documentation/2019/reference-manual/algorithms/expanded-ensemble.html
    # The use of the alpha is not completly right from a theretical point of view for this method
    assert len(temperatures) == len(energies), "temperatures and energies must have the same number of elements"
    # Deffining the coefficient matrix
    M = 1*np.eye(len(temperatures), k = 1, dtype=int)
    np.fill_diagonal(M, -1)
    M[-1] = np.ones(len(temperatures), dtype = int)
    # Solution Vector
    b = np.zeros(len(temperatures))
    # Initialize b[0] dependign if alpha was selected
    if alpha:
        gamma = alpha/(1-alpha)
        b[0] = 0.5*(tools.beta(temperatures[1]) - tools.beta(temperatures[0]))*(energies[0] + energies[1]) - np.log(gamma*(len(temperatures) - 1))
    else:
        b[0] = 0.5*(tools.beta(temperatures[1]) - tools.beta(temperatures[0]))*(energies[0] + energies[1])

    for i in range(1, len(energies)):
        try:
            b[i] = 0.5*(tools.beta(temperatures[i + 1]) - tools.beta(temperatures[i]))*(energies[i] + energies[i + 1])
        # I reached the last element of the list so I need to add cero to the array but I already initialized on zero.
        except:
            break
    # Solve the system of equations
    # Calculate the weights:
    weights = np.linalg.solve(M, b)
     
    # Convert if needed to GROMACS weights
    if weight_GROMACS_format:
        minimum = weights.min()
        if minimum > 0:
            weights -= minimum
        else:
            weights += np.abs(minimum)
    return weights

def get_energy(edr, annealing_times, energy_type = 'Potential', out_fig = 'energy_distribution.svg'): # Could be Total-Energy
    """The two input is a list that could be obtained from mdp.annealing3
    The structure is the following
    annealing_times = [t1, t2, t3, t4, t5, t6]
    These mean:
    From mdp.annealing, we also obtain the list
    annealing_temperatures = [T1, T1, T2, T2, T3, T3]
    Tis two list are connected and mean:
    from time t1 to t2 the temperature was constant, therefore from here we should get the average energy for T1.
    Then from t2 to t3 the temperature increase from T1 to T2
    Then from t3 to t4 the temperature is constant at T2,therefore from here we should get the average energy for T2.
    And so on and so for.

    Args:
        annealing_temperatures (_type_): _description_
        annealing_times (_type_): _description_
    """
    fig, ax = plt.subplots(figsize = (16,9))
    data = pd.DataFrame()
    xvg_tmp_file = tempfile.NamedTemporaryFile(suffix='.xvg')
    energy = []
    iterator = range(0, len(annealing_times)-1, 2)

    for state, index in tqdm.tqdm(enumerate(iterator), total=len(iterator)):#enumerate(iterator):# # the calculation is per pair of times, beetween the first to time the temperature was keep constant, then the system was heated and repeated again.
        run = tools.run(f"export GMX_MAXBACKUP=-1; echo {energy_type} | gmx energy -f {edr} -b {annealing_times[index]} -e {annealing_times[index + 1]} -o {xvg_tmp_file.name} | grep \'{energy_type.replace('-',' ')}\'")
        energy.append(float(run.stdout.split()[-5]))
        """
        Energy                      Average   Err.Est.       RMSD  Tot-Drift
        -------------------------------------------------------------------------------
        Potential                -1.30028e+06         --     1682.1   -2422.24  (kJ/mol)
        Total Energy                -952595         --    2606.81    -3688.3  (kJ/mol)
        """
        # Getting the histograms and checking for the same len in all intervals
        if state == 0:
            data[state] = xvg.XVG(xvg_tmp_file.name).data[:,1]
        else:
            xvg_data = xvg.XVG(xvg_tmp_file.name).data[:,1]
            if xvg_data.shape[0] > data.shape[0]:
                data[state] = xvg_data[:data.shape[0]]
            else:
                data = data.iloc[:xvg_data.shape[0]]
                data[state] = xvg_data


    print(data)
    sns.histplot(data = data, element='poly', stat = 'probability', axes = ax)
    ax.set(
        xlabel = f'{energy_type} [kJ/mol]',
        ylabel = 'Probability',
        title = f'Distribution of {energy_type}')
    # plt.show()
    fig.savefig(out_fig)
    return energy

def get_weights_from_log(log, plot = False):
    """This function read a log file obtained from simulated tempering and returns a numpy array
    of shape (number of time points, number of states, 2)
    The first axis of the array refears to each output on the log file.
    The second axis refears to the temperature states
    And the last one to the counts and value of the corresponded weights
    E.g., if we would like to know the value of the weight on the second output of the log file for the temperature state 5.
    >>> from palmiche.simutemp import weights
    >>> time, info_weights = weights.get_weights_from_log('tempering.log')
    >>> info_weights[1,4,1]
    387.69299

    Get the final weights:
    >>> info_weights[-1,:,1]
    array([   0.     ,   95.78638,  194.94769,  291.73962,  387.69299,
        484.79767,  563.36194,  609.72424,  656.0224 ,  702.54724,
        747.32648,  793.23291,  838.1333 ,  880.62061,  922.65308,
        966.77673, 1008.16711, 1051.17847, 1091.65161, 1133.02808,
       1175.10168, 1212.64014, 1249.82324, 1289.07568, 1329.84229,
       1366.30273])

    Args:
        log (str): path to the log file

    Returns:
        numpy.array: The weights information.
    """
    with open(log, 'r') as f:
        log_file = f.readlines()

    i = 0
    time = []
    weights_info = []
    weights_0 = []
    while i < len(log_file):
        if 'init-lambda-weights[' in log_file[i]:
            weights_0.append(float(log_file[i].split('=')[-1]))

        if 'MC-lambda information' in log_file[i]:
            # Finding the time
            for j in range(i,0,-1):
                if log_file[j].startswith('           Step           Time'):
                    j += 1
                    time.append(float(log_file[j].split()[-1]))
                    break
            # Finding the weight
            weights_info_tmp = []
            i += 3
            while log_file[i] != '\n':
                split = log_file[i].split()
                count = int(split[2])
                weight = float(split[3])
                weights_info_tmp.append((count, weight))
                i += 1
            weights_info.append(weights_info_tmp)
        i += 1
    # Add weights at t = 0, because the counts are all 0 and I delate the entrances with total count 0 in next lines,
    # What i could do is put 1 in the initial temperature
    time.insert(0,0)
    weights_info.insert(0,list(zip([1] + (len(weights_0) - 1)*[0], weights_0)))

    #Converting to array
    time = np.array(time)
    weights_info = np.array(weights_info)
    # Some times (I don't know why) GROMACS reset all the weights and all the counts are 0. We need to eliminate those points
    sum_of_weights = weights_info[:,:,0].sum(axis = 1)
    time = time[sum_of_weights != 0]
    weights_info = weights_info[sum_of_weights != 0]
    sum_of_weights = sum_of_weights[sum_of_weights != 0]


    if plot:
        dir = os.path.dirname(log)
        fig, axes = plt.subplots(2, figsize = (16,9), sharex=True)
        NUM_COLORS = weights_info.shape[1]
        cm = plt.get_cmap('viridis')#gist_rainbow viridis
        for axe in axes:
            axe.set_prop_cycle('color', [cm(1.*j/NUM_COLORS) for j in range(NUM_COLORS)])

        probability = weights_info[:,:,0] / sum_of_weights[:,np.newaxis]
        for j in range(weights_info.shape[1]):
            #axes[0].plot(time, weights_info[:,j,0], label = str(j))
            axes[0].plot(time, probability[:,j], label = str(j))
            axes[1].plot(time, weights_info[:,j,1])
        
        fig.legend(loc = 'lower center', ncol = int(weights_info.shape[1] / 2))
        axes[0].set(
            xlim = (time.min(), time.max()),
            ylim = (0,1),
            ylabel = 'Probability',
            )
        axes[1].set(
            xlabel = 'Time [ps]',
            ylabel = 'Weight values'
        )
        #plt.show()
        fig.savefig(os.path.join(dir,'weights_progression.svg'), bbox_inches="tight")

        # Plotting the violin plot of the weights
        df = pd.DataFrame()
        for j in range(weights_info.shape[1]):
            #df[temperatures[j]] = weights_info[:,j,1]
            df[j] = weights_info[:,j,1]
        # Set up the matplotlib figure
        sns.set_theme(style="whitegrid")
        fig, ax = plt.subplots(figsize=(25, 25))

        # Draw a violinplot with a narrower bandwidth than the default
        sns.violinplot(data=df, palette="Set3", bw=.2, cut=1, linewidth=1)
        # The plot is not over the actual temperatures, the temperatures ara only labels
        ax.plot(range(len(weights_info[0,:,1])), weights_info[0,:,1], '-o',  label = 'Initial weights')
        ax.set(
            title = 'Weights per state over the entire simulation',
            xlabel = 'Sate',
            ylabel = 'Weight',
        )
        plt.legend()
        sns.despine(left=True, bottom=True)
        #plt.show()
        fig.savefig(os.path.join(dir,'weights_per_state.svg'), bbox_inches="tight")
        sns.reset_defaults()

    return time, weights_info


if __name__ == '__main__':
    print("You are in the submodule weights of palmiche.simutemp")
    # a = get_weights_from_log('/home/ale/mnt/smaug/TEST/500_H2O/tempering.log')
    # print(a)

    # o = get_energy('/home/ale/mnt/smaug/TEST/500_H2O/tempering.edr', [0,20,25,30,50,200, 230,240,563,600,700,800,900,950])
    # print(o)

    # a = PandeWeights([100,400,500], [20,17,256], alpha = 0.0001)
    # print(a)
    # temperatures = list(np.linspace(300, 501, 400))
    # energies = range(len(temperatures))
    # target_temperatures = [300, 325, 350, 375, 400, 425, 450, 475, 500]
    # a = AleWeights(temperatures, energies, target_temperatures, alpha = 2)
    # print(a)
    # path = '/home/ale/mnt/smaug/TEST/water/AleWeights'#'/home/ale/mnt/smaug/TEST/500_H2O/'#'/home/ale/mnt/smaug/TEST/water/AleWeights'
    # temp_list = range(300, 332, 2)
    # temperature = xvg.XVG(os.path.join(path,'temperature.xvg')).data
    # total_energy = xvg.XVG(os.path.join(path,'total_energy.xvg')).data
    # # mask = (temperature[:,1] >= min(temp_list) - 2) & (temperature[:,1] <= max(temp_list) + 2)
    # # temperature = np.array(temperature[mask])
    # # total_energy = np.array(total_energy[mask])
    # print(temperature[:,1])
    # A = AleWeights(
    #     temperature[:,1],
    #     total_energy[:,1],
    #     temp_list,
    #     alpha = None,
    #     weight_GROMACS_format = True,
    #     integrator = 'trapezoid')
    # print(A)
    _, w = get_weights_from_log('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_P/7e27/sys_MMV007839_Cell_891_SP_param/windows/00085/production.log', plot = True) # AleWeights