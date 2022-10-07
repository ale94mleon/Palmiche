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
from argparse import ArgumentParser
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
        


if __name__=='__main__':
    #Getting the arguments
    #setupUmbrella.py --disfile(-d) distance_file --interval(-i) sample_interval --template(-t) template_file
    parser = ArgumentParser(description='Getting umbrella windows (This script assumes an increases in the COM)')
    parser.add_argument('-d', '--disfile', dest='disfile',
                    help="""The name of distance file in the format n*2, where
                    the two columns are the the integer identify of the 
                    conformer and the center of mass""")
    parser.add_argument('-i1', '--interval1', dest='interval1',
                    help="""
                    The main interval distance between each windows. 
                    """)
    parser.add_argument('-i2', '--interval2', dest='interval2',
                    help="""
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
                    dismissed.
                    """)
    parser.add_argument('-n2', '--number2', dest='number2', default=1,
                    help="""
                    Control of how many windows are created with the 
                    secondary interval close to the window of lower COM. 
                    Default is 1.
                    """)
    parser.add_argument('-t', '--template', dest='template',
                    help="""
                    The template for the execution of GROMACS
                    """)
    parser.add_argument('-r', '--resume', dest='resume', default="resume.dat",
                    help="""
                    Resume of the windows distance
                    """)
    args = parser.parse_args()
    
    for item in args.__dict__:
        if item not in ['interval2','template'] :
            if args.__dict__[item] == None:
                exit(f"The argument \"{item}\" is missing. The right way to call the script is (in \"[]\" the optional arguments): {__usage__}")
    
    
    # Figure out which frames to use
    distance_table = readDistanceFile(args.disfile)
    #Becasue first = True, the first conformation is always added
    sampled_indexes = sampleDistances(distance_table,float(args.interval1),first=True)
    
    #If we want to add extra windows near to the window with the lowest COM:
    #    interval2 controls the width and number2 how many.
    if args.interval2:
        sampled_indexes += sampleDistances(distance_table,float(args.interval2))[1:int(args.number2)+1]
        
    sampled_indexes = list(sorted(set(sampled_indexes))) #to remove the duplicate items
    
    # If any template files were specified, use them to make frame-specific
    # output
    if args.template != None:
        print("Creating frame-specific output for {args.template}:")
        for i in sampled_indexes:
            frame = distance_table[i][0]
            createOutputFile(args.template,frame,search_string="XXX")    
    
    # Print out summary of the frames we identified 
    to_print = sorted([distance_table[index] for index in sampled_indexes], key=lambda x:x[1])

    out = ["%10s%10s%10s\n" % ("frame","dist","d_dist")]#El Conformer_ID, reaction coordinate distance, diference with the former distance
    for i, item in enumerate(to_print):
        frame = item[0]
        dist = item[1]
        if i == 0:
            delta_dist = "%10s" % "NA" 
        else:
            prev_dist = to_print[i-1][1]
            delta_dist = "%10.3f" % (dist - prev_dist)
        
        out.append("%10i%10.3f%s\n" % (frame,dist,delta_dist))
    
    print("".join(out))
    with open(args.resume, 'w') as o:
        o.writelines(out)
