#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Tue Jun 15 16:07:59 2021
Author        : Alejandro Martínez León
Mail          : [alejandro.martinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""
__description__ = \
"""
Take a file containing frame number vs. distance, then identify frames that 
sample those distances at approximately sample_interval.  Optionally takes an
arbitrary number of template files.  The program searches the contents of these
files for a search string (by default XXX), replaces the search string with the
frame number, then writes out each file with a unique, frame-specific filename.
This last feature means that this script can be used to automatically generate
input files to run each umbrella sampling simulation. 

v. 0.1: In this version I (Alejadnro) made modification in the way to select
the umbrella windows (explain why) and also some other small changes respect
to the version 0.1 of Michael J. Harms [harmsm@gmail.com].
"""
__usage__ = "setupUmbrella.py --disfile(-d) distance_file --interval1(-i1) main_sample_interval [--interval2(-i2) secondary_sample_interval --number(-n2) number_of_secondary_interval --template(-t) template_file]"

import numpy as np
import pandas as pd
from tools import backoff

 
def readDistanceFile(distance_file):
    with open(distance_file, 'r') as f:
        lines = f.readlines()
    return [(int(line.split()[0]), float(line.split()[1])) for line in lines]

def sampleDistances(table_distance,sample_interval, first = False):
    """
    Go through the distances list and sample frames at sample_interval.  
    Appropriate samples are identified by looking forward through the
    distances to find the one that is closest to current_distance +
    sample_interval.  
    """
    
    distance = [d[1] for d in table_distance]

    windows = []
    windows_index = []
    miss_center = []
    for center in np.arange(min(distance), max(distance) + sample_interval, sample_interval):
        onward = [abs(com-center) for com in distance]
        min_index = onward.index(min(onward))
        if center - sample_interval/2 <= distance[min_index]\
        and distance[min_index] <= center + sample_interval/2\
        and distance[min_index] not in windows:
            windows.append(distance[min_index])
            windows_index.append(min_index)
        else:
            miss_center.append(center)
    
           
    if miss_center: print(f"Warning: For the interval {sample_interval} the following center of windows were not found: {miss_center}")
    if first:
        if 0 not in windows_index:
            windows_index.insert(0,0)
    return windows_index

def createOutputFile(template_file,frame_number,search_string="XXX"):
    """
    Look for instances of the search string in a template file and replace with
    the frame number in a new file.
    """

    out_file = f"frame-{frame_number}_{template_file}"

    # Prompt the user before wiping out an existing file
    backoff(out_file)
    
    # Read the contents of the template file
    with open(template_file,'r') as f:
        file_contents = f.read()

    # Write out the template file contents, replacing all instances of 
    # the search string with the frame number
    with open(out_file,'w') as f:
        f.write(file_contents.replace(search_string,str(frame_number)))
        


def windows_creator(disfile, interval1, interval2 = 0, num_windows_interval2 = 1, template = None):
    """
    

    Parameters
    ----------
    disfile : TYPE: string
        DESCRIPTION.
        The name of distance file in the format n*2, where
        the two columns are the the integer identify of the 
        conformer and the center of mass
    
    interval1 : TYPE: float
        DESCRIPTION.
        The main interval distance between each window. In nm
    
    interval2 : TYPE: float
        DESCRIPTION. The default is 0. It will not used
        This interval is used in case do you want to add some other
        windows respect to the first window (the one with lower COM, 
        which is not mandatory the starting configuration for the 
        pulling, could be another one. For example, if you want 
        ,close to the first window, one window of 0.1 and your
        main interval is 0.2 (controlled by interval1). Then you 
        need to set --interval2(-i2) 0.1. If you wnat more than one
        windows close to the first one you need to change 
        number2, that by default is 1. If the new windows is 
        already included with the main interval, then it will be 
        dismissed    

    num_windows_interval2 : TYPE, optional: float
        DESCRIPTION. The default is 1.
        How many windows close to the first with the second interval as distance
   
    template : TYPE, optional: str
        DESCRIPTION. The default is None.
        The template for the execution of GROMACS


    Returns
    -------
    A DataFrame with the resume of the selected windows with three columns:
        ["Conf_ID","Distance","Delta_Distance"]

    """   
    # Figure out which frames to use
    distance_table = readDistanceFile(disfile)
    #Becasue first = True, the first conformation is always added
    sampled_indexes = sampleDistances(distance_table,float(interval1),first=True)
    
    #If we want to add extra windows near to the window with the lowest COM:
    #    interval2 controls the width and number2 how many.
    if interval2:
        sampled_indexes += sampleDistances(distance_table,float(interval2))[1:int(num_windows_interval2)+1]
        
    sampled_indexes = list(sorted(set(sampled_indexes))) #to remove the duplicate items
    
    # If any template files were specified, use them to make frame-specific
    # output
    if template != None:
        print("Creating frame-specific output for {template}:")
        for i in sampled_indexes:
            frame = distance_table[i][0]
            createOutputFile(template,frame,search_string="XXX")    
    
    # Print out summary of the frames we identified 
    to_print = sorted([distance_table[index] for index in sampled_indexes], key=lambda x:x[1])

    out = []
    for i, item in enumerate(to_print):
        frame = item[0]
        dist = item[1]
        if i == 0:
            delta_dist = "NaN" 
        else:
            prev_dist = to_print[i-1][1]
            delta_dist = (dist - prev_dist)
        
        out.append([frame,dist,delta_dist])
    
    to_return = pd.DataFrame(data=out, columns = ["Conf_ID","Distance","Delta__Distance(i)-Distance(i-1)"])
    
    return to_return
