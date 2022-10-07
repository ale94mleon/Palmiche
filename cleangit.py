#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Thu Oct 1 00:18:59 2021
Author        : Alejandro Martínez León.
Mail          : [alejandro.martinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""
from palmiche.tools import tools
from glob import glob
import os
for (root, dirs, files) in os.walk("."):
    for file in files:
        if file.endswith(".pyc"):
            pyc_path = os.path.join(root, file)
            #print(pyc_path)
            tools.run(f"git rm -f {pyc_path}")
[tools.run(f"git rm -f {item} -r") for item in ["./examples", "__pycache__", "./umbrella/TEST"]]

print("hello")