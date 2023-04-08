#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from palmiche.utils import tools, xvg
import argparse 


class POT:
    """
    This class is used to calculate the flat bottom potential of a trajectory where this kind of potential was used. See:
    https://manual.gromacs.org/documentation/2019/reference-manual/functions/restraints.html
    """
    def __init__(self, typepot, coords,**keywords):#force
        self.typepot = typepot
        #self.force = force
        self.coords = coords
        self.keywords = keywords
        POTENTIALS = {"flat_bottom": np.vectorize(self.flat_bottom)}
        self.pot = POTENTIALS[self.typepot](self.coords)
    def flat_bottom(self, array):
        if array < self.keywords['ref']:
            return 0.0
        else:
            return 0.5*self.keywords['k']*(array - self.keywords['ref'])**2

        #self.pot = [-0.5*(item[0] - self.keywords['ref'])*item[1] for item in zip(self.coord, self.force)]

    def plot(self, time):
        plt.plot(time, self.pot)
        plt.title(f"Potential.")
        plt.xlabel("Time")
        plt.ylabel("V(t)")
        plt.show()
      
       
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='A test program.')
    parser.add_argument("coord", help = "The index of the coord (column index in the cvg file).", type = int)
    # parser.add_argument("-f", 
    #                     dest="pullf",
    #                     default = "pull_pullf.xvg",
    #                     help = "The pullf file. Default is pull_pullf.xvg",
    #                     type = str)
    parser.add_argument("-x", 
                        dest="pullx",
                        default = "pull_pullx.xvg",
                        help = "The pullx file. Default is pull_pullx.xvg",
                        type = str)   
    parser.add_argument("-r", 
                        dest="ref",
                        default = 0,
                        help = "The reference value for the potential. Default is 0",
                        type = float)
    parser.add_argument("-k", 
                        dest="k",
                        default = 0,
                        help = "The reference value for the potential. Default is 0",
                        type = float)
    parser.add_argument("-p", 
                        dest="pot",
                        default = "flat_bottom",
                        help = "The potential name. Default is flat_bottom. Could be selected: flat_bottom, ...",
                        type = str)  
    args = parser.parse_args()

    #f = xvg.XVG(args.pullf).data[:,args.coord]
    x = xvg.XVG(args.pullx).data[:,args.coord]
    t = xvg.XVG(args.pullx).data[:,0]
    k = POT(args.pot, x, ref = args.ref, k = args.k).plot(t)
