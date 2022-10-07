#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Wed Sep 15 23:16:58 2021
Author        : Alejandro Martínez León
Mail          : [alejandro.martinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""
import os
from glob import glob
from palmiche.utils import tools, heavyh, mdp
import matplotlib.pyplot as plt

# class MDP:
#     def __init__(self, *args, dt = 0.004, time = 250, xtc_numb_frame = 2000):
        
#         """
        
    
#         Parameters
#         ----------
#         *args : TYPE, tup
#             DESCRIPTION. Any possible tuple of two members: (key, value)
#         dt : TYPE, optional, float
#             DESCRIPTION. The default is 0.004. The interval of integration
#         time : TYPE, optional, float
#             DESCRIPTION. The default is 250. The integration time in ns
#         xtc_numb_frame : TYPE, optional. int
#             DESCRIPTION. The default is 2000. How many frmaes would you like to have in the trajectory file
        
#         dt and time are used to calculate the numebr of step
#         and xtc_numb_frame to calculate the steps elapsed.
#         If you add some argument, the calculated value will change by the provided ones.

    
#         """
#         self.args = args
#         self.dt = dt
#         self.time = time
#         self.xt_numb_frame =xtc_numb_frame
#         self.nsteps = int(time * 1000 / dt)
#         self.steps_elapsed = int(self.nsteps / xtc_numb_frame) 
#         self.template = f"""integrator              = md
# dt                      = {self.dt}
# nsteps                  = {self.nsteps}		
# nstxout                 = 0				
# nstvout                 = 0				
# nstfout                 = 0				
# nstcalcenergy           = {int(self.steps_elapsed / 2)}
# nstenergy               = {self.steps_elapsed}		
# nstlog                  = {self.steps_elapsed}				
# nstxout_compressed      = {self.steps_elapsed}				
# cutoff_scheme           = Verlet
# nstlist                 = 20
# rlist                   = 1.0				
# vdwtype                 = Cut-off
# vdw_modifier            = Potential-shift-Verlet		
# rvdw_switch             = 0				
# rvdw                    = 1.0				
# coulombtype             = pme
# rcoulomb                = 1.0				
# epsilon_r               = 1
# epsilon_rf              = 1
# tcoupl                  = v-rescale
# tc_grps                 = SOLU MEMB SOLV
# tau_t                   = 1.0 1.0 1.0  
# ref_t                   = 303.15 303.15 303.15
# pcoupl                  = c-rescale
# pcoupltype              = semiisotropic
# tau_p                   = 5.0
# compressibility         = 4.5e-5  4.5e-5
# ref_p                   = 1.0     1.0
# constraints             = h-bonds
# constraint_algorithm    = LINCS
# continuation            = yes
# nstcomm                 = 100
# comm_mode               = linear
# comm_grps               = SOLU MEMB SOLV
# refcoord_scaling        = com"""
#         self.string = ""
#         self.dict_template = {}
#         self.parse()
#     def parse(self):
        
#         for line in self.template.split("\n"):
#             if not line.startswith(';'): 
#                 split = [x.strip() for x in line.split("= ")]
#                 try:
#                     self.dict_template[split[0]] = split[1]
#                 except:
#                     pass
            
#         for option in self.args:
#             print(f"Added: {option[0]} = {option[1]}\n")
#             self.dict_template[option[0].replace("-", "_")] = option[1]
        
#         for key in self.dict_template:
#             self.string += f"{key:<40} = {self.dict_template[key]:<15}\n"
        
#     def write(self, name = 'production.mdp'):
#         with open(name, "w") as out:
#             out.write(self.string)



def mass_factor_analysis(*replica_paths, log = "md.log", info_write = True, info_graph = True):
    summary = {}
    for replica in replica_paths:
        for root, dirs, files in os.walk(replica):
            if log in files:
                mass_factor = float(os.path.basename(root).replace("_", "."))
                check = tools.CHECK(os.path.join(root, log))
                
                if mass_factor not in summary:
                    if check.performance:
                        summary[mass_factor] = {"ended_simulations":1, "lower_performance":check.performance, "higher_performance":check.performance}
                    else:
                        summary[mass_factor] = {"ended_simulations":0, "lower_performance":0,  "higher_performance":0}                   

                else:
                    if check.performance:
                        summary[mass_factor]["ended_simulations"] += 1
                        if check.performance < summary[mass_factor]["lower_performance"]:
                            summary[mass_factor]["lower_performance"] = check.performance
                        if check.performance > summary[mass_factor]["higher_performance"]:
                            summary[mass_factor]["higher_performance"] = check.performance
    
    if info_write:
        with open("info_HMR_mass_factor.txt", "w") as info:
            info.write(f"{'mass_factor':<15}{'lower_performance':<20}{'higher_performance':<20}{'ended_simulations':<20}\n")
            for mass_factor in sorted(summary.keys()):
                info.write(f"{mass_factor:<15}{summary[mass_factor]['lower_performance']:<20}{summary[mass_factor]['higher_performance']:<20}{summary[mass_factor]['ended_simulations']:<20}\n")
    
    if info_graph:
        width = 0.50
        x_labels = []
        y_axis = []
        bar_label = []
        for mass_factor in sorted(summary.keys()):
            x_labels.append(str(mass_factor))
            y_axis.append(summary[mass_factor]['ended_simulations'])
            bar_label.append(f"l:{int(summary[mass_factor]['lower_performance'])}-h:{int(summary[mass_factor]['higher_performance'])}")
        
        fig, ax = plt.subplots(figsize=(15, 5))
        p1 = ax.bar(x_labels, y_axis, width)
        ax.set_yticks(range(1,len(replica_paths)+1))
        ax.set_ylabel('Ended simulations')
        ax.set_xlabel('H-mass factors')
        ax.set_title('H-mass repartioning')
        #ax.bar_label(p1, labels = bar_label, rotation='vertical', label_type='center')
        ax.bar_label(p1, labels = bar_label)
        plt.show()
        fig.savefig('myfig.svg')
    return summary


        
    


def dict_gen(path = "./assemble_mini_equi", FF_protein = "amber99sb-star-ildn.ff",
             FF_lipid = "Slipids_2020.ff" , FF_ligand = "GAFF2.ff", topol = "topol.top",
             toppar = "toppar", conf = "step6.6_equilibration.gro", ndx = "index.ndx"):
    """
    THos fucntion demand to have a tree shape:
        path/receptor1/ligand1 and so on
        path/receptor2/ligand456 and so on and with the fiels or directories
        difained in the function inside the last directorie (ligand

    Parameters
    ----------
    path : TYPE, optional
        DESCRIPTION. The default is "./assemble_mini_equi".
    FF_protein : TYPE, optional
        DESCRIPTION. The default is "amber99sb-star-ildn.ff".
    FF_lipid : TYPE, optional
        DESCRIPTION. The default is "Slipids_2020.ff".
    FF_ligand : TYPE, optional
        DESCRIPTION. The default is "GAFF2.ff".
    topol : TYPE, optional
        DESCRIPTION. The default is "topol.top".
    toppar : TYPE, optional
        DESCRIPTION. The default is "toppar".
    conf : TYPE, optional
        DESCRIPTION. The default is "step6.6_equilibration.gro".

    Returns
    -------
    results : TYPE
        DESCRIPTION.

    """
    path = os.path.abspath(path) 
    dict_paths = {}
    for receptor in tools.list_if_dir(path):
        for ligand in tools.list_if_dir(os.path.join(path, receptor)): #In this way if the receptor where teasted to different ligands  (een so, they must be in the deffinition of the ligands_dict) this will tak it intop account
            dict_paths.setdefault(receptor,{})[ligand] = {"FF_protein":os.path.join(path, receptor, ligand, FF_protein),
                                                       "FF_lipid":os.path.join(path, receptor, ligand, FF_lipid),
                                                       "FF_ligand":os.path.join(path, receptor, ligand, FF_ligand),
                                                       "topol":os.path.join(path, receptor, ligand, topol),
                                                       "toppar":os.path.join(path, receptor, ligand, toppar),
                                                       "conf":os.path.join(path, receptor, ligand, conf),
                                                       "ndx":os.path.join(path, receptor, ligand, ndx),
                                                       }
    return dict_paths

def HMR(dict_paths, mdp_args, path_out = "HMR_3_0", heavyh_massfactor = 3.0,
        heavyh_revert = 0, heavyh_suff = "_heavyH"):
    
    mdp_obj = mdp.MDP(**mdp_args)
    
    cwd = os.getcwd()
    path_out = os.path.abspath(path_out)
    transfer_path = os.path.abspath("transfer")
    tools.makedirs(transfer_path)
    os.chdir(transfer_path)
    for receptor in dict_paths:
        for ligand in dict_paths[receptor]:
            dest_path = os.path.join(path_out, receptor, ligand)
            dest_toppar_path = os.path.join(path_out, receptor, ligand, "toppar")
            tools.makedirs(dest_path)
            tools.makedirs(dest_toppar_path)
            
            [tools.cp(item, transfer_path) for item in [os.path.join(dict_paths[receptor][ligand]["toppar"], "*"), dict_paths[receptor][ligand]["topol"]]]
            
            topol_name = os.path.basename(dict_paths[receptor][ligand]["topol"])
            
            with open(topol_name, "r") as f:
                data = f.read()
            data = data.replace("toppar/","")
            with open(topol_name, "w") as f:
                f.write(data)
            #print(data)
            heavyh.run_heavyh(topol = topol_name, massfactor = heavyh_massfactor,
                              revert = heavyh_revert, suff = heavyh_suff)
            #Getting the new itp created and moving to the corresponfing directory
            heavy_itp = glob(f"*{heavyh_suff}.itp") 
            [tools.mv(item, dest_toppar_path) for item in heavy_itp]
            
            topol_name = f"{topol_name.split('.')[0]}{heavyh_suff}.top"
            with open(topol_name, "r") as f:
                data = f.read()
            #!!These repleacment is system depended
            data = data.replace("#include \"topol_Protein","#include \"toppar/topol_Protein")
            data = data.replace("#include \"LI","#include \"toppar/LI")
            data = data.replace("#include \"POPC","#include \"toppar/POPC")
            with open(topol_name, "w") as f:
                f.write(data)
            
            #include "topol_Protein
            #include "LI
            
            tools.mv(topol_name, dest_path)
            
            #If we eliminate the heavyh_suff from the name we will get the original name
            heavy_itp = [file.replace(heavyh_suff, "") for file in heavy_itp]
            rest_of_itp = list(set(glob("*.itp")) - set(heavy_itp))
            [tools.mv(item, dest_toppar_path) for item in rest_of_itp]
            
            #Cleaning tansfer
            tools.rm(os.path.join(transfer_path, "*"))
            
            [tools.cp(item, dest_path, r = True) for item in [dict_paths[receptor][ligand]["FF_protein"],
                                                   dict_paths[receptor][ligand]["FF_lipid"],
                                                   dict_paths[receptor][ligand]["FF_ligand"],
                                                   dict_paths[receptor][ligand]["conf"],
                                                   dict_paths[receptor][ligand]["ndx"]
                                                   ]]
            mdp_obj.write(os.path.join(dest_path, "production.mdp"))
                     
    tools.rm(transfer_path, r = True)
    os.chdir(cwd)
    
if __name__ == '__main__':
    print("We are in HMR_maker")