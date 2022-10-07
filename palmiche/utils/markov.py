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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

class MARKOV:
    def __init__(self, realization):
        self.realization = realization
        #>> matrix = Markov(realization)
        #>> matrix.compute_transition_matrix() will update the matrix if a different realization is passed
        self.compute_transition_matrix()
    
    def compute_transition_matrix(self):
        #Getting the states check what method to use for sorted the states
        if any(isinstance(n, str) for n in self.realization):
            states = sorted(set(self.realization), key=lambda x:str(x))
        else:
            states = sorted(set(self.realization))

        #Translate states to numerical indexes
        transdict = dict(zip(states, range(len(states))))

        # Converting the realization to numbers
        realization_to_work = [transdict[item] for item in self.realization]
        
        # Creating the matrix
        self.transition_matrix = np.zeros([len(states)]*2)
        for (i,j) in zip(realization_to_work,realization_to_work[1:]):
            self.transition_matrix [i][j] += 1
        
        # Convert to probability
        for row in self.transition_matrix :
            s = row.sum()
            if s > 0:
                row /= s

        # Adding labels
        self.transition_matrix = pd.DataFrame(self.transition_matrix, columns=states, index = states)
    
    def plot(self, out = 'transition_matrix.svg', occupancy = True):
        sns.set_theme(style="whitegrid")
        if occupancy:
            fig, ax = plt.subplots(nrows=1, ncols=2, figsize = (25,10))
            sns.heatmap(self.transition_matrix, ax = ax[0], vmin=.0, vmax=.5, center=.3, annot=True, cbar_kws={'label': 'Transition Probability'})
            df = pd.DataFrame(Counter(self.realization), index = [0])
            df /= df.values.sum()
            # Rearrange the columns
            cols = df.columns.tolist()
            if any(isinstance(n, str) for n in cols):
                cols = sorted(cols, key=lambda x:str(x))
            else:
                cols = sorted(cols)
            sns.barplot(data = df[cols], ax = ax[1])
            ax[0].set(
                xlabel = 'States',
                ylabel = 'States')
            ax[1].set(
                xlabel = 'States',
                ylabel = 'Probability')
        else:
            fig, ax = plt.subplots(figsize = (12,10))
            sns.heatmap(self.transition_matrix, ax = ax, vmin=.0, vmax=.5, center=.3, annot=True, cbar_kws={'label': 'Transition Probability'})
            ax.set(
                xlabel = 'States',
                ylabel = 'States'
                )
        fig.savefig(out, bbox_inches="tight")
        sns.reset_defaults()

        
# def transition_matrix(realization):
#     #The convertion is in case different types of data are passed
#     realization = np.array(realization)
    
#     #Getting the states
#     states = np.sort(np.unique(realization))

#     #Translate states to numerical indexes
#     transdict = dict(zip(states, range(len(states))))

#     # Converting the realization to numbers
#     realization_to_work = [transdict[item] for item in realization]
    
#     # Creating the matrix
#     M = np.zeros([states.shape[0]]*2)
#     for (i,j) in zip(realization_to_work,realization_to_work[1:]):
#         M[i][j] += 1
    
#     # Convert to probability
#     for row in M:
#         s = row.sum()
#         if s > 0:
#             row /= s
    
#     # Adding labels
#     M = pd.DataFrame(M, columns=states, index = states)
#     return M

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import seaborn as sns
    from palmiche.utils import xvg
    realization = xvg.XVG('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_P/7e27/sys_MMV007839_Cell_891_SP_param/windows/00001/production.xvg').data[:,1]
    # realization = [1,1,1,2,2,'B',4,4,8,8,55,47,'A', 'B', 'A']
    # realization = [1,1,1,2,2,4,4,8,8,55,47,5,5,1,1,2,3,10,12,13,13,13,13,15,15,45,47,48,47,45,46,49,48,47]
    plt.plot(realization)
    plt.show()
    A = MARKOV(realization.astype(int))
    A.plot('p.pdf')
