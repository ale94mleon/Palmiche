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
import copy, warnings
import numpy as np


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

class MDP:
    """This is a wrap around the mdp file of GROMACS, the units are the used by GROMACS:
        time: ps
        distance: nm
        etc...
    """
    def __init__(self, type = 'production', **user_keywords) -> None:
        self.type = type
        self.user_keywords = user_keywords
        self.default_keywords = {
            'minimization':{
                'define':'',
                'integrator':'steep',
                'emtol':1000.0,
                'nsteps':5000,
                'nstlist':10 ,
                'cutoff_scheme':'Verlet',
                'rlist':1.0,
                'vdwtype':'Cut-off',
                'vdw_modifier':'Potential-shift-Verlet',# Force-switch
                'rvdw_switch':0,
                'rvdw':1.0,
                'coulombtype':'pme',
                'rcoulomb':1.0,
                'epsilon_r':1,
                'epsilon_rf':1,
                'constraints':'h-bonds',
                'constraint_algorithm':'LINCS',
            },
            'equilibration':{
                # Examples of define
                # define = -DPOSRES -DPOSRES_FC_BB=1000.0 -DPOSRES_FC_SC=100.0 -DPOSRES_FC_LIPID=1000.0 -DDIHRES -DDIHRES_FC=1000.0 -DPOSRES_LIG=0.0 -DFLAT_BOTOM_k=600.0 -DFLAT_BOTOM_r=0.5
                'define':'',
                'integrator': 'md',
                'dt': 0.002,
                'tinit':0,
                'nsteps':None,
                'nstcomm':100,
                #Output parameters
                'nstxout':0,
                'nstvout':0,
                'nstfout':0,
                'nstcalcenergy':100,
                'nstenergy':5000,
                'nstlog':None,
                'nstxout_compressed':None,
                #Single-range cutoff scheme
                'cutoff_scheme':'Verlet',
                'nstlist':50,
                'rlist':1.0	,
                'vdwtype':'Cut-off',
                'vdw_modifier':'Potential-shift-Verlet',
                'rvdw_switch':0,
                'rvdw':1.0,
                'coulombtype':'pme',
                'rcoulomb':1.0,
                'epsilon_r':1,
                'epsilon_rf':1,
                #Temperature coupling
                'Tcoupl':'v-rescale',
                'tc_grps':'System',
                'tau_t':1.0,
                'ref_t':303.15,
                #Pressure coupling
                'pcoupl':'c-rescale',#Used for both, equilibration and production
                # semiisotropic is usefull for membrane simulation, and the project associated with this script was on this direction
                # But if you have, for example protein in water
                # pcoupltype=isotropic,tau_p=5.0,compressibility=4.5e-5,ref_p=1.0 could be an option
                'pcoupltype':'semiisotropic',
                'tau_p':5.0,
                'compressibility':'4.5e-5  4.5e-5',
                'ref_p':'1.0  1.0',
                #Bond parameters
                'constraint_algorithm':'lincs',
                'constraints':'h-bonds',
                'continuation':'yes',
                #Generate velocities is off
                'gen_vel':'no',
                #Periodic boundary conditions are on in all directions
                'pbc':'xyz',
                #Long-range dispersion correction
                'DispCorr':'EnerPres',
                'refcoord_scaling':'com',
            },
            'production':{
                # Examples of define
                # define = -DPOSRES -DPOSRES_FC_BB=1000.0 -DPOSRES_FC_SC=100.0 -DPOSRES_FC_LIPID=1000.0 -DDIHRES -DDIHRES_FC=1000.0 -DPOSRES_LIG=0.0 -DFLAT_BOTOM_k=600.0 -DFLAT_BOTOM_r=0.5
                'define':'',
                'integrator': 'md',
                'dt': 0.002,
                'tinit':0,
                'nsteps':None,
                'nstcomm':100,
                #Output parameters
                'nstxout':0,
                'nstvout':0,
                'nstfout':0,
                'nstcalcenergy':100,
                'nstenergy':5000,
                'nstlog':None,
                'nstxout_compressed':None,
                #Single-range cutoff scheme
                'cutoff_scheme':'Verlet',
                'nstlist':50,
                'rlist':1.0	,
                'vdwtype':'Cut-off',
                'vdw_modifier':'Potential-shift-Verlet',
                'rvdw_switch':0,
                'rvdw':1.0,
                'coulombtype':'pme',
                'rcoulomb':1.0,
                'epsilon_r':1,
                'epsilon_rf':1,
                #Temperature coupling
                'Tcoupl':'v-rescale',
                'tc_grps':'System',
                'tau_t':1.0,
                'ref_t':303.15,
                #Pressure coupling
                'pcoupl':'c-rescale',#Used for both, equilibration and production
                # semiisotropic is usefull for membrane simulation, and the project associated with this script was on this direction
                # But if you have, for example protein in water
                # pcoupltype=isotropic,tau_p=5.0,compressibility=4.5e-5,ref_p=1.0 could be an option
                'pcoupltype':'semiisotropic',
                'tau_p':5.0,
                'compressibility':'4.5e-5  4.5e-5',
                'ref_p':'1.0  1.0',
                #Bond parameters
                'constraint_algorithm':'lincs',
                'constraints':'h-bonds',
                'continuation':'yes',
                #Generate velocities is off
                'gen_vel':'no',
                #Periodic boundary conditions are on in all directions
                'pbc':'xyz',
                #Long-range dispersion correction
                'DispCorr':'EnerPres',
                'refcoord_scaling':'com',
            },
        }
        # This condition holds only for the case that equilibration and production is used, and only on those scenarios must be provided.
        # Time is in [ns]
        # This section is for automation on some variables that depends on the time, dt, and desired frames of simulation
        # Checking for the keywords time and xtc_numb_frame in case that equilibration or production
        if self.type in ['equilibration','production']:
            if 'dt' in self.user_keywords:
                dt = self.user_keywords['dt']
            else:
                dt = self.default_keywords[self.type]['dt']

            if 'time' not in self.user_keywords:
                warnings.warn(f"Because {self.type} was selected and you did not provided the user_keyword 'time' (simulation time in ps); "\
                    f"The value 1000 ps will be used.")
                self.user_keywords['time'] = 1000

            if 'xtc_numb_frame' not in self.user_keywords:
                warnings.warn(f"Because {self.type} was selected and you did not provided the user_keyword 'xtc_numb_frame' (Number of frames of the trajectory); "\
                    f"The value 10 will be used. ")
                self.user_keywords['xtc_numb_frame'] = 10
            # In case that nstcalenergy is provided, it will be used to calculate the frequency of the output of the log file (nstlog)
            if 'nstcalcenergy' in self.user_keywords:
                self.default_keywords[self.type]['nstcalcenergy'] = self.user_keywords['nstcalcenergy']

            nsteps = int(self.user_keywords['time'] / dt)
            steps_elapsed = int(nsteps /  self.user_keywords['xtc_numb_frame'])

            # Deffine nstlog, in case that there are not to many steps (small simulations usually for testing) set to nstcalcenergy
            nstlog = min([int(nsteps / 2000), steps_elapsed]) # Here we get at least 2000 points. 
            if nstlog <= self.default_keywords[self.type]['nstcalcenergy']:
                nstlog = self.default_keywords[self.type]['nstcalcenergy']

            self.default_keywords[self.type]['dt'] = dt
            self.default_keywords[self.type]['nsteps'] = nsteps
            self.default_keywords[self.type]['nstlog'] =  nstlog
            self.default_keywords[self.type]['nstxout_compressed'] = steps_elapsed

        self.keywords = {}
        self.keywords['General'] = copy.deepcopy(self.default_keywords[self.type])
        
        for key in self.user_keywords:
            # Do not take into account these two keys that are not mdp valid options
            if key not in  ['time', 'xtc_numb_frame']:
                # Convert all '-' to '_' in order that the key is identified on self.keywords
                key = key.replace('-', '_')
                self.keywords['General'][key] = self.user_keywords[key]
                if key not in self.default_keywords:
                    print(f'Added keyword: {key} = {self.user_keywords[key]}')
        


    def pull(self, dist2pull, ligands2pull, flat_bottom_init = 0.4, flat_bottom_k = 500, orient_restrain_init = 45, orient_restrain_k = 8274.41, orient_restrain_vec = (0,0,-1), **pullkeywords):
        """
        This will genrate the pull section with cylinder flat bottom restrain and orientation restrain
        orient_restrain_init = 0.349 # 20 degrees 20/180*3.141
        orient_restrain_k = 1308.36 # in Kj/mol*rad^2, leads to sigma = 2.5 degrees

        Parameters
        ----------
        dist2pull : TYPE, optional, float
            DESCRIPTION. How long we want to pull in nm.

        ligands2pull : TYPE, optional, list of strings
            DESCRIPTION. The list of names used in the index file for the ligands to pull
            if you just want to use the function to modify mdp options
        *args : TYPE, tup
            DESCRIPTION. Any possible tuple of two members: (key, value)

        Returns
        -------
        None.

        """
        ncoords = 3*len(ligands2pull)
        ngroups = 4*len(ligands2pull)
        pull_code = {
            # Pull code
            'pull':'yes',
            'pull_ncoords':ncoords,# reaction coordinates
            'pull_ngroups':ngroups,# different groups
            'pull_nstxout':50,# for the production 125
            'pull_nstfout':50,# for the production 125
            # REACT_COORD_DEF
            'pull_group1_name':'LIA_CLOSE_AA',# These are the residues at some nm from the ligand (should have been defined in automation) at 1.2 nm
            'pull_group2_name':'LIA',# This is the ligand
            'pull_coord1_type':'umbrella',# harmonic potential
            'pull_coord1_geometry':'direction',
            'pull_coord1_dim':'N N Y',
            'pull_coord1_vec':'0.0 0.0 -1.0',
            'pull_coord1_groups':'1 2',
            'pull_coord1_start':'no',#define initial COM distance > 0,
            'pull_coord1_init':0,# see for production
            'pull_coord1_k':2000,# kJ mol^-1 nm^-2
        }
        if self.type == 'minimization':
            pull_code['pull_coord1_rate'] = 0.0
        else:
            pull_code['pull_coord1_rate'] = dist2pull / (self.user_keywords['time'])
        #Addding and changing default options
        for key in pullkeywords:
            key = key.replace('-', '_')
            if key not in pull_code:
                print(f'Added keyword: {key} = {pull_code[key]}')
            pull_code[key] = pullkeywords[key]

        #Spliting the dict into two dictionary, the first part and the reaction coordinates definition
        #pull2build it will be used for construct this section
        pull2build = {}
        pull_coords_section = {}
        for key in pull_code:
            if key.startswith('pull_group') or key.startswith('pull_coord'):
                pull_coords_section[key] = pull_code[key]
            else:
                pull2build[key] = pull_code[key]

        for (i,ligand) in enumerate(ligands2pull):

            for key in pull_coords_section:
                if key == 'pull_group1_name':
                    pull2build[f'pull_group{2*i+1}_name'] = f"{ligand}_CLOSE_AA"
                elif key == 'pull_group2_name':
                    pull2build[f'pull_group{2*i+2}_name'] = ligand
                elif key == 'pull_coord1_groups':
                    pull2build[f'pull_coord{i+1}_groups'] = f"{2*i+1} {2*i+2}"
                else:
                    pull2build[key.replace('coord1', f'coord{i+1}')] = pull_coords_section[key]

        self.keywords['Pull Code'] = dict()
        for key in pull2build:
            self.keywords['Pull Code'][key] = pull2build[key]

        # Here I will define the cylinder flat bottom potential
        #lateral cylindrical flat-bottomed restraint, so the molecule does not miss the pore entrance\n"
        pull_flat_bottom_section = {
            'pull_coord1_type':'flat-bottom',
            'pull_coord1_geometry':'distance',
            'pull_coord1_dim':'Y Y N',
            'pull_coord1_groups':'1 2',
            'pull_coord1_start':'no',
            'pull_coord1_init':flat_bottom_init,
            'pull_coord1_rate': 0.0,
            'pull_coord1_k':flat_bottom_k,
        }
        pull2build_flat_bottom = dict()
        for (i,ligand) in enumerate(ligands2pull):

            for key in pull_flat_bottom_section:
                if key == 'pull_coord1_groups':
                    pull2build_flat_bottom[f'pull_coord{i + 1 + len(ligands2pull)}_groups'] = f"{2*i + 1} {2*i + 2}"
                else:
                    pull2build_flat_bottom[key.replace('coord1', f'coord{i + 1 + len(ligands2pull)}')] = pull_flat_bottom_section[key]
        # Here because I am also creating the string, but a more pythonic way to update the dict is just
        # self.keywords.update(pull2build_flat_bottom)
        self.keywords['Lateral Cylindrical Flat-bottomed Restraint'] = dict()
        for key in pull2build_flat_bottom:
            self.keywords['Lateral Cylindrical Flat-bottomed Restraint'][key] = pull2build_flat_bottom[key]

        # Here I will define the orientation restrains
        orientation_section = {
            'pull_group1_name':'LIA_GROUP1',
            'pull_group2_name':'LIA_GROUP2',
            'pull_coord1_type':'flat-bottom',
            'pull_coord1_geometry':'angle-axis',
            'pull_coord1_dim':'Y Y Y',
            'pull_coord1_vec':' '.join([str(xi) for xi in orient_restrain_vec]),
            'pull_coord1_groups':'1 2',
            'pull_coord1_start':'no',
            'pull_coord1_init':orient_restrain_init,
            'pull_coord1_rate': 0.0,
            'pull_coord1_k':orient_restrain_k,
        }

        pull2build_orientation_restrain = {}
        for (i,ligand) in enumerate(ligands2pull):

            for key in orientation_section:
                if key == 'pull_group1_name':
                    pull2build_orientation_restrain[f'pull_group{2*i + 1 + 2*len(ligands2pull)}_name'] = f"{ligand}_GROUP1"
                elif key == 'pull_group2_name':
                    pull2build_orientation_restrain[f'pull_group{2*i + 2 + 2*len(ligands2pull)}_name'] = f"{ligand}_GROUP2"
                elif key == 'pull_coord1_groups':
                    pull2build_orientation_restrain[f'pull_coord{i + 1 + 2*len(ligands2pull)}_groups'] = f"{2*i + 1 + 2*len(ligands2pull)} {2*i + 2 + 2*len(ligands2pull)}"
                else:
                    pull2build_orientation_restrain[key.replace('coord1', f'coord{i + 1 + 2*len(ligands2pull)}')] = orientation_section[key]

        self.keywords['Orientation Restraint'] = dict()
        for key in pull2build_orientation_restrain:
            self.keywords['Orientation Restraint'][key] = pull2build_orientation_restrain[key]


    def annealing1(self, NumbGroups, temp, heat_fraction = 0.25, nstenergy_same_as_nstcalcenergy = False):
        """
        This is usefull for AleWeights
        Generate the annealing section of the MDP.
        The annealing strategy is the folloging:
        Heat the system during heat_fraction*self.time at the begining and aftter cold the system till
        the end of the simulation. For all the groups will be used 'single' as type of annealing and will be
        used the same temperatures.
        E.g:
            NumGroups = 2
            temp = [200,210,220]
            self.time = 25000 ps
            heat_fraction = 0.25

        OUT:
            annealing = 'single single'
            annealing_npoints = '5 5'
            annealing_time = '0 3125 6250 15625 25000 0 3125 6250 15625 25000'
            annealing_temp = '200 210 220 210 200 200 210 220 210 200'


        Args:
            NumbGroups (int): The number of temperature couplead and to be used in the Annealing
            temp (list): List of temperature to use the
            heat_fraction (float, optional): How much of the simulation time will be spent in the heating process. Defaults to 0.25.
        """
        # In case that asked, set nstenergy to nstcalcenergy value
        if nstenergy_same_as_nstcalcenergy:
            self.keywords['General']['nstenergy'] = self.keywords['General']['nstcalcenergy']
        temp = np.array(temp)
        annealing ={
            'annealing':' '.join(NumbGroups*['single']),
            'annealing_npoints':' '.join(NumbGroups*[str(2*len(temp) - 1)]), # Two times temp because heat and cold and minus 1 because the cold has one lest temperature
        }

        annealing_time =  np.hstack([
            np.linspace(0, heat_fraction*self.user_keywords['time'], len(temp)),
            np.linspace(heat_fraction*self.user_keywords['time'], self.user_keywords['time'], len(temp))[1:]
        ]).round(2)
        annealing_time = ' '.join([str(t) for t in annealing_time])
        annealing['annealing_time'] = ' '.join([annealing_time for i in range(NumbGroups)])

        annealing_temp = np.hstack([
            temp,
            np.flip(temp[:-1])
        ]).round(2)
        annealing_temp = ' '.join([str(t) for t in annealing_temp])
        annealing['annealing_temp'] = ' '.join([annealing_temp for i in range(NumbGroups)])

        self.keywords['Simulated Annealing'] = dict()
        for key in annealing:
            self.keywords['Simulated Annealing'][key] = annealing[key]


    def annealing2(self, GroupBool, temp_list, heating_frac = 0.25, cooling_frac = 0.5):
        """
        In this case we are able to specify which groups will change the temperature.
        Generate the annealing section of the MDP.
        The annealing strategy is the folloging:
        Heat the system during heat_fraction*self.time at the begining and aftter cold the system till cooling_frac*self.time. Then mantain the target temperature. This is useful for equilibration.
        the end of the simulation. For all the groups will be used 'single' as type of annealing and will be
        used the same temperatures.

        @@@@@@@ Mejorar el ejemplo@@@@@@@
        E.g:
            NumGroups = 2
            temp = [200,210,220]
            self.time = 25000 ps
            heat_fraction = 0.25

        OUT:
            annealing = 'single single'
            annealing_npoints = '5 5'
            annealing_time = '0 3125 6250 15625 25000 0 3125 6250 15625 25000'
            annealing_temp = '200 210 220 210 200 200 210 220 210 200'


        Args:
            NumbGroups (int): The number of temperature couplead and to be used in the Annealing
            temp (list): List of temperature to use the
            heat_fraction (float, optional): How much of the simulation time will be spent in the heating process. Defaults to 0.25.
        """
        temp_list = np.array(temp_list)
        annealing ={
            'annealing':' '.join(len(GroupBool)*['single']),
            'annealing_npoints':' '.join(len(GroupBool)*[str(2*len(temp_list) - 1)]), # Two times temp_list because heating and cooling and minus 1 because the cooling has one lest temperature
        }

        annealing_time =  np.hstack([
            np.linspace(0, heating_frac*self.user_keywords['time'], len(temp_list)),
            np.linspace(heating_frac*self.user_keywords['time'], (heating_frac + cooling_frac)*self.user_keywords['time'], len(temp_list))[1:]
        ]).round(2)
        annealing_time = ' '.join([str(t) for t in annealing_time])
        annealing['annealing_time'] = ' '.join([annealing_time for i in GroupBool])

        annealing_temp = np.hstack([
            temp_list,
            np.flip(temp_list[:-1])
        ]).round(2)
        # In case that we dont need to dynamically change the temeprature of one group we only fix to tyhe target value that is the first temerature
        annealing_temp_static = ' '.join([str(annealing_temp[0]) for t in annealing_temp])
        # For the groups that we want to change the temperature and generate a possible better sampling
        annealing_temp_dynamic = ' '.join([str(t) for t in annealing_temp])
        annealing['annealing_temp'] = ' '
        for test in GroupBool:
            if test:
                 annealing['annealing_temp'] += annealing_temp_dynamic + ' '
            else:
                annealing['annealing_temp'] += annealing_temp_static + ' '

        self.keywords['Simulated Annealing'] = dict()
        for key in annealing:
            self.keywords['Simulated Annealing'][key] = annealing[key]

    def annealing3(self, NumbOfGroups, temp, heating_time = 100, constant_temp_time = 2000):
        """
        Generate the annealing section of the MDP.
        The annealing strategy is the following:

        T3             _____
                       |
        T2        ____|
                 |       
        T1  ____|
        This is used as input of PandeWeights method

        Heat the system during heat_fraction*self.time at the begining and aftter cold the system till
        the end of the simulation. For all the groups will be used 'single' as type of annealing and will be
        used the same temperatures.
        @@@@Mejorar el ejemplo, cambiar la variable temp por temperatures
        E.g:
            NumGroups = 2
            temp = [200,210,220]
            self.time = 25000 ps
            heat_fraction = 0.25

        OUT:
            annealing = 'single single'
            annealing_npoints = '5 5'
            annealing_time = '0 3125 6250 15625 25000 0 3125 6250 15625 25000'
            annealing_temp = '200 210 220 210 200 200 210 220 210 200'


        Args:
            NumbGroups (int): The number of temperature couplead and to be used in the Annealing
            temp (list): List of temperature to use the
            heat_fraction (float, optional): How much of the simulation time will be spent in the heating process. Defaults to 0.25.
        """
        temp = np.array(temp)
        annealing ={
            'annealing':' '.join(NumbOfGroups*['single']),
            'annealing_npoints':' '.join(NumbOfGroups*[str(2*len(temp))]), # Two times temp because for each temperature there are tow points where the temperature stays constant..
        }
        annealing_time = [0]
        for _ in temp: # Do not add to the last point a heaing time
            annealing_time.append(annealing_time[-1] + constant_temp_time)
            annealing_time.append(annealing_time[-1] + heating_time)
        # Deleting the last heating time
        annealing_time.pop()
        # Changing the simulaiton time
        # Getting the dt
        if self.type in ['equilibration','production']:
            if 'dt' in self.user_keywords:
                dt = self.user_keywords['dt']
            else:
                dt = self.default_keywords[self.type]['dt']
        else:
            raise ValueError(f"type = {self.type} and must be production or equilibration for anneling simulation")
        nsteps = int(annealing_time[-1]/ dt)
        steps_elapsed = int(nsteps /  self.user_keywords['xtc_numb_frame'])

        nstlog = min([int(nsteps / 2000), steps_elapsed]) # Here we get at least 2000 points. 
        if nstlog <= self.keywords['General']['nstcalcenergy']: # In case that the output frequancy is higher than the calculation of the energy.
            nstlog = self.keywords['General']['nstcalcenergy']

        self.keywords['General']['dt'] = dt
        self.keywords['General']['nsteps'] = nsteps
        # self.keywords['General']['nstcalcenergy'] = 100 # int(steps_elapsed / 2) # I need to think in a more cleaver way to do this. But the main problem is that here I need a good statistic of the enrgy, so, I need to calculate and export more frequentlly.
        self.keywords['General']['nstenergy'] = 10 * self.keywords['General']['nstcalcenergy']
        self.keywords['General']['nstlog'] = nstlog
        self.keywords['General']['nstxout_compressed'] = steps_elapsed
        
        for key in self.user_keywords:
            # Do not take into account these two keys that are not mdp valid options
            if key not in  ['time', 'xtc_numb_frame']:
                # Convert all '-' to '_' in order that the key is identified on self.keywords
                key = key.replace('-', '_')
                self.keywords['General'][key] = self.user_keywords[key]


        # here we 
        annealing_time_str = ' '.join([str(t) for t in annealing_time])
        annealing['annealing_time'] = ' '.join([annealing_time_str for _ in range(NumbOfGroups)])
       
        annealing_temp = np.repeat(temp,2)
        annealing_temp_str = ' '.join([str(t) for t in annealing_temp])
        annealing['annealing_temp'] = ' '.join([annealing_temp_str for _ in range(NumbOfGroups)])

        self.keywords['Simulated Annealing'] = dict()
        for key in annealing:
            self.keywords['Simulated Annealing'][key] = annealing[key]
        return annealing_temp, annealing_time


    def simulated_tempering(self, temp_list, init_lambda_weights, **st_keywords):

        # Getting the lambda values for simulated tempering on GROMACS from a list of temperatures states.
        # In this case simulated_tempering_scaling must be equal to linear.
        temp_list = np.array(temp_list)
        T_MAX = temp_list.max()
        T_MIN = temp_list.min()
        temperature_lambdas = (temp_list - T_MIN) / (T_MAX - T_MIN)

        # Change the integration method
        self.keywords['General']['integrator'] = 'md-vv'        
        # # Increase the frequency of the energy calculation and energy output
        # self.keywords['General']['nstcalcenergy'] = int(self.keywords['General']['nstcalcenergy'] / 2)
        # self.keywords['General']['nstenergy'] = 2 * self.keywords['General']['nstcalcenergy']
               
        # Control  the frequency of nstdhl output.
        nstdhdl = self.keywords['General']['nstcalcenergy']
        if 'Pull Code' in self.keywords:
            nstdhdl = self.keywords['Pull Code']['pull_nstxout']

        tempering = {
            # Free energy variables
            'free_energy': 'expanded',
            'init_lambda_state': 0,
            'nstdhdl': nstdhdl,
            'temperature_lambdas':' '.join([str(tl) for tl in temperature_lambdas]),
            'init_lambda_weights':' '.join([str(il) for il in init_lambda_weights]),
            # Simulated tempering variables
            'simulated_tempering':'yes',
            'simulated_tempering_scaling':'linear',
            'sim_temp_low':T_MIN,
            'sim_temp_high':T_MAX,
            # Expanded ensemble variables
            'nstexpanded':self.keywords['General']['nstcalcenergy'],
            'lmc_stats':'wang-landau',
            'lmc_move':'metropolis', #lmc-mc-move = metropolis-transition
            'lmc_weights_equil':'no',
            # 'weight_equil_number_all_lambda':-1,
            # 'weight_equil_number_samples':-1,
            # 'weight_equil_number_steps':-1,
            # 'weight_equil_wl_delta':-1,
            # 'weight_equil_count_ratio':-1,
            # Seed for Monte Carlo in lambda space
            'wl_scale':0.999999,
        }
        #Addding and changing default options
        for key in st_keywords:
            key = key.replace('-', '_')
            # This is a safe check in order to do not change this parameter, becasue the lambda were calculated following the inear strategy.
            if key == 'simulated_tempering_scaling' and st_keywords != 'linear':
                raise ValueError(f"You provided for '{key}' a value of '{st_keywords[key]}'; however only 'linear' is accepted.")
            if key not in tempering:
                print(f'Added keyword: {key} = {tempering[key]}')
            tempering[key] = st_keywords[key]

        self.keywords['Simulated Tempering'] = dict()
        for key in tempering:
            self.keywords['Simulated Tempering'][key] = tempering[key]

    def string(self):
        mdp_str = ''
        for first_key in self.keywords:
            if first_key != 'General':
                mdp_str += f"\n; {first_key}\n"

            for second_key in self.keywords[first_key]:
                mdp_str += f"{second_key:<40} = {self.keywords[first_key][second_key]}\n"
        return mdp_str


    def write(self, name = 'mdp.mdp'):
        with open(name, "w") as out:
            out.write(self.string())

if __name__ == '__main__':
    import json
    tc_grps = 'SOLU MEMB SOLV LIA LIB LIC LID LIE'
    mdp = MDP(
        type = 'equilibration',
        dt = 0.004,
        time = 100,
        xtc_numb_frame = 2000,
        define = '-DPOSRES -DPOSRES_FC_BB=1000.0 -DPOSRES_FC_SC=100.0 -DPOSRES_FC_LIPID=1000.0 -DDIHRES -DDIHRES_FC=1000.0 -DPOSRES_LIG=0.0',
        tc_grps = 'SOLU MEMB SOLV LIA LIB LIC LID LIE',
        tau_t = '1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0',
        ref_t = '303.15 303.15 303.15 303.15 303.15 303.15 303.15 303.15'
        )
    #mdp.pull(4.5,['LI'+c for c in 'ABCDE'])
    #mdp.annealing1(1,np.arange(303.15,355,2),heat_fraction=0.25)
    tc_grps = np.array(tc_grps.split())
    #mask = (tc_grps != 'MEMB') & (tc_grps != 'SOLV')
    #mdp.annealing2(mask,np.arange(303.15,303.15+6,2), heating_frac=0.25, cooling_frac=0.5)
    a, b = mdp.annealing3(2, np.arange(303.15,303.15+6,2))
    mdp.pull(4.5, [4,5,6,7,8])
    #mdp.simulated_tempering([300,54,587], [6,5,8])
    print(mdp.string())
    # print(a,b)
    # print(mdp.keywords['Pull Code']['pull_nstfout'])
    # print(mdp.keywords['Simulated Tempering']['nstexpanded'])
    # mdp.keywords['Simulated Tempering']['nstexpanded'] = mdp.keywords['Pull Code']['pull_nstfout']
    # mdp.keywords['Simulated Tempering']['nstdhdl'] = mdp.keywords['Pull Code']['pull_nstfout']
    # print(mdp.keywords['Pull Code']['pull_nstfout'])
    # print(mdp.keywords['Simulated Tempering']['nstexpanded'])
    # print(mdp.string())


