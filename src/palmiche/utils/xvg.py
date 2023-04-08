#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shutil, os
import re

def backoff(file_path):
    """
    

    Parameters
    ----------
    file : TYPE string
        DESCRIPTION: The name or the path f     or the specific file

    Returns
    -------
    None.
    If the file already exist. it will made a back up to ./#{file}.{str(i)}#,
    Where i is an integer.
    """
    basname = os.path.basename(file_path)
    dirname = os.path.dirname(file_path)
    if os.path.exists(file_path):
        new_basname = basname
        i = 1
        while(os.path.exists(os.path.join(dirname, new_basname))):
            new_basname = f"./#{basname}.{str(i)}#"
            i += 1
        print(f"Back Off! I just backed up {file_path} to {os.path.join(dirname, new_basname)}")
        shutil.copy2(file_path, os.path.join(dirname, new_basname))

class XVG:
    def __init__(self, file):
        self.file = file
        self.title = ""
        self.type = ""
        self.legend = []
        self.xlable = ""
        self.ylable = ""
        self.data = []
        self.multidata = []
        self.parser()

    def parser(self):
        with open(self.file, "r") as f:
            for line in f.readlines():
                if line.startswith("#"):
                    continue
                if line.startswith("@"):
                    if " title" in line:
                        self.title = line.split("title")[1].strip().replace("\"", "")
                    elif "TYPE" in line:
                        self.type = line.split()[1]
                    elif "xaxis  label" in line:
                        self.xlabel = line.split("xaxis  label")[1].strip().replace("\"", "")
                    elif "yaxis  label" in line:
                        self.ylabel = line.split("yaxis  label")[1].strip().replace("\"", "")
                    # Looking for the labels
                    item = re.search('@ s[0-9]+ legend', line) # Match @ s25 legend
                    if item:
                        self.legend.append(re.search("\"([^;]*)\"", line).group(1)) # Match the name of the label
                    continue
                if line.startswith("&"):
                    self.data = np.asarray(self.data, dtype=float)
                    self.multidata.append(self.data)
                    self.data = []
                    continue

                self.data.append(line.split())
        self.data = np.asarray(self.data, dtype=float)

    def write(self, out_file = None):
        if out_file:
            pd.DataFrame(self.data).to_string(out_file, index = False, header = False, float_format= '%.5f')
        else:
            backoff(self.file)
            pd.DataFrame(self.data).to_string(self.file, index = False, header = False, float_format= '%.5f')

    def quickgraph(self):
        fig, ax = plt.subplots()
        NUM_COLORS = self.data.shape[1]
        cm = plt.get_cmap('gist_rainbow')
        ax.set_prop_cycle('color', [cm(1.*i/NUM_COLORS) for i in range(NUM_COLORS)])
        for i in range(1, self.data.shape[1]):
            ax.plot(self.data[:,0], self.data[:,i], label = (self.legend[i-1]))
        ax.set(
            title = self.title,
            xlabel = self.xlabel,
            ylabel = self.ylabel,
        )
        fig.legend()
        plt.show()
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
                        help = "xvg file",
                        dest='xvg',
                        type=str)
    args = parser.parse_args()
    xvg = XVG(args.xvg)
    #print(xvg.xlabel, xvg.title, xvg.ylabel)
    xvg.quickgraph()

                    