#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__description__ = \
"""
Take a file containing frame number vs. distance, then identify frames that 
sample those distances at approximately window_width.  Optionally takes an
arbitrary number of template files.  The program searches the contents of these
files for a search string (by default XXX), replaces the search string with the
frame number, then writes out each file with a unique, frame-specific filename.
This last feature means that this script can be used to automatically generate
input files to run each umbrella sampling simulation. 

v. 0.1: In this version I (Alejandro) made modifications in the way to select
the umbrella windows (below explained why) and also some other small changes respect
to the first version of Michael J. Harms [harmsm@gmail.com].
"""
__usage__ = "setupUmbrella.py --disfile(-d) distance_file --interval1(-i1) main_window_width [--interval2(-i2) secondary_window_width --number(-n2) number_of_secondary_interval --template(-t) template_file]"

import os, numpy as np
from argparse import ArgumentParser
 
def readDistanceFile(distance_file):
    with open(distance_file, 'r') as f:
        lines = f.readlines()
    return [(int(line.split()[0]), float(line.split()[1])) for line in lines]

def sampleDistances(table_distance, window_width, first = False, **keywords):
    """
    Go through the distances list and sample frames at window_width.  
    Appropriate samples are identified by looking forward through the
    distances to find the one that is closest to current_distance +
    window_width.  


    Args:
        table_distance (list): This list is returned by the function readDistanceFile()
        window_width (float): This is the width of each window.py
        first (bool, optional): In case that the first frame want to be for sure included. Defaults to False.
        **keywords (optional): minimum and maximum. 
            In case you want to specified the range for the umbrella, if not is specified minimum and/or maximum. The real
            extreme value(s) of table_distance will be used.
            In the case of maximum, if you want to pick the value you need to provide, for example:
                maximum + window_width, because if not the np.arange will not take into account the last number.py. If you don't specify
                the maximum value the function will add automatically window_width to maximum.

    Returns:
        [list]: [List of window indexes]
    """
    
    distance = [d[1] for d in table_distance]
    #Here we are taken the interval to be sampled
    try:
        minimum = keywords["minimum"]
    except:
        minimum = min(distance)
    try:
        maximum = keywords["maximum"]
    except:
        maximum = max(distance) + window_width


    windows = []
    windows_index = []
    miss_center = []
    for center in np.arange(minimum, maximum, window_width):
        onward = [abs(com-center) for com in distance]
        min_index = onward.index(min(onward))
        if center - window_width/2 <= distance[min_index]\
        and distance[min_index] <= center + window_width/2\
        and distance[min_index] not in windows:
            #print(center, distance[min_index])
            windows.append(distance[min_index])
            windows_index.append(min_index)
        else:
            miss_center.append(center)
    
           
    if miss_center:
        print(f"Warning: For the interval {window_width} the following center of windows were not found: {miss_center}")
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
    if os.path.exists(out_file):
        answer = input(f"{out_file} exists!  Overwrite (y|n)?")
        answer = answer[0].lower()
        if answer != "y":
            return None
    
    # Read the contents of the template file
    with open(template_file,'r') as f:
        file_contents = f.read()

    # Write out the template file contents, replacing all instances of 
    # the search string with the frame number
    with open(out_file,'w') as f:
        f.write(file_contents.replace(search_string,str(frame_number)))

def main(disfile, window_width, first = False, template = None, txt_summary = 'conf_windows_summary.txt',  **keywords):
    # Figure out which frames to use
    distance_table = readDistanceFile(disfile)
    #Because first = True, the first conformation is always added
    sampled_indexes = sampleDistances(distance_table,float(window_width),first=first, **keywords)

    # If any template files were specified, use them to make frame-specific
    # output
    if template != None:
        print(f"Creating frame-specific output for {template}:")
        for i in sampled_indexes:
            frame = distance_table[i][0]
            createOutputFile(template,frame,search_string="XXX")    
    
    # Print out summary of the frames we identified 
    to_print = sorted([distance_table[index] for index in sampled_indexes], key=lambda x:x[1])

    
    out_string = f"{'window_ID':<15}{'frame':<10}{'dist':<10}{'delta_dist':<10}\n"
    windows = {}
    for i, item in enumerate(to_print):
        window_ID = i+1
        frame = item[0]
        dist = item[1]
        if i == 0:
            delta_dist = f"{'NaN':<10}" 
        else:
            prev_dist = to_print[i-1][1]
            delta_dist = f"{(dist - prev_dist):<10.3}"
        
        out_string += f"{window_ID:<15}{frame:<10}{dist:<10.3}{delta_dist}\n"
        windows[window_ID] = {'frame': frame, 'dist': dist, 'delta_dist': delta_dist.strip()}
    
    if txt_summary:
        with open(txt_summary, 'w') as o:
            o.write(out_string)
    return windows

     


def mainxxxx(disfile, interval1, interval2 = 0, num_windows_interval2 = 1, template = None):
    """
    This is only the main function, the description is forward. type:
        -h to get help

    """    
    
    # Figure out which frames to use
    distance_table = readDistanceFile(disfile)
    #Because first = True, the first conformation is always added
    sampled_indexes = sampleDistances(distance_table,float(interval1),first=True)
    
    #If we want to add extra windows near to the window with the lowest COM:
    #    interval2 controls the width and num_windows_interval2 how many.
    if interval2:
        sampled_indexes += sampleDistances(distance_table,float(interval2))[1:int(num_windows_interval2)+1]
        
    sampled_indexes = list(set(sampled_indexes)) #to remove the duplicate items
    
    # If any template files were specified, use them to make frame-specific
    # output
    if template != None:
        print("Creating frame-specific output for {template}:")
        for i in sampled_indexes:
            frame = distance_table[i][0]
            createOutputFile(template,frame,search_string="XXX")    
    
    # Print out summary of the frames we identified 
    to_print = sorted([distance_table[index] for index in sampled_indexes], key=lambda x:x[1])

    
    out_string = f"{'window_ID':<15}{'frame':<10}{'dist':<10}{'delta_dist':<10}\n"
    windows = {}
    for i, item in enumerate(to_print):
        window_ID = i+1
        frame = item[0]
        dist = item[1]
        if i == 0:
            delta_dist = f"{'NA':<10}" 
        else:
            prev_dist = to_print[i-1][1]
            delta_dist = f"{(dist - prev_dist):<10.3}"
        
        out_string += f"{window_ID:<15}{frame:<10}{dist:<10.3}{delta_dist}\n"
        windows[window_ID] = {'frame': frame, 'dist': dist, 'delta_dist': delta_dist.strip()}
    
    with open('resume_conf_windows.txt', 'w') as o:
        o.write(out_string)
    return windows

if __name__=='__main__':
    #Getting the arguments
    #setupUmbrella.py --disfile(-d) distance_file --interval(-i) window_width --template(-t) [template_file1 template_file2 ...]
    parser = ArgumentParser(description='Getting umbrella windows (This script assumes an increases in the COM)')
    parser.add_argument('-d', '--disfile', dest='disfile',
                    help="""The name of distance file in the format n*2, where
                    the two columns are the the integer identify of the 
                    conformer and the center of mass""")
    parser.add_argument('-i1', '--interval1', dest='interval1',
                    help="""
                    The main interval distance between each windows. In nm 
                    """)
    parser.add_argument('-i2', '--interval2', dest='interval2', default = 0,
                    help="""
                    The default is 0,  not usedT.his interval is used in case
                    do you want to add some other windows respect to the first
                    window (the one with lower COM, which is not mandatory the
                    starting configuration for the pulling, could be another one.
                    For example, if you want, close to the first window, one 
                    window of 0.1 and your main interval is 0.2 (controlled by interval1).
                    Then you need to set --interval2(-i2) 0.1. If you want more than one
                    windows close to the first one you need to change 
                    num_windows_interval2, that by default is 1. If the new windows is 
                    already included with the main interval, then it will be 
                    dismissed.
                    """)
    parser.add_argument('-n2', '--num_windows_interval2', dest='num_windows_interval2', default=1,
                    help="""
                    Control of how many windows are created with the 
                    secondary interval close to the window of lower COM. 
                    Default is 1.
                    """)
    parser.add_argument('-t', '--template', dest='template',
                    help="""
                    The template for the execution of GROMACS (a shell script)
                    """)
    args = parser.parse_args()
    
    for item in args.__dict__:
        if item not in ['interval2','template']:
            if args.__dict__[item] == None:
                raise SyntaxError(f"The argument \"{item}\" is missing. The right way to call the script is (in \"[]\" the optional arguments): {__usage__}")
    
    main(args.disfile, args.interval1, interval2 = args.interval2, num_windows_interval2 = args.num_windows_interval2, template = args.template)
