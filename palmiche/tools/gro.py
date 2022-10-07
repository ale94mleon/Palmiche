#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Wed Jun  9 15:32:10 2021
Author        : Alejandro Martínez León
Mail          : [alejandro.martinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""
from palmiche.tools import tools
import numpy as np

class InputERROR(Exception):
    pass

class Atom:
    def __init__(self, line):
        #The last is not take
        self.resNumber = int(line[:5]) #residue number (5 positions, integer)
        self.resName = line[5:10].strip() #residue name (5 characters)
        self.atomName = line[10:15].strip() #atom name (5 characters)
        self.atomNumber = int(line[15:20]) #atom number (5 positions, integer)
        #position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
        self.x = float(line[20:28]) 
        self.y = float(line[28:36])
        self.z = float(line[36:44])
        #velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)
        try: 
            self.vX = float(line[44:52])
        except:
            self.vX = 0.0000
        try: 
            self.vY = float(line[52:60])
        except:
            self.vY = 0.0000   
        try: 
            self.vZ = float(line[60:68])
        except:
            self.vZ = 0.0000
        
    def __getitem__(self, key):
        return self.__dict__[key]

class VECTOR():
    def __init__(self,line):
        split = line.split()
        try: 
            self.v1x = float(split[0])
        except: 
            self.v1x = 0.0
        try: 
            self.v2y = float(split[1])
        except: 
            self.v2y = 0.0 
        try: 
            self.v3z = float(split[2])
        except: 
            self.v3z = 0.0        
        try: 
            self.v1y = float(split[3])
        except: 
            self.v1y = 0.0        
        try: 
            self.v1z = float(split[4])
        except: 
            self.v1z = 0.0  
        try: 
            self.v2x = float(split[5])
        except: 
            self.v2x = 0.0      
        try: 
            self.v2z = float(split[6])
        except: 
            self.v2z = 0.0 
        try: 
            self.v3x = float(split[7])
        except: 
            self.v3x = 0.0  
        try: 
            self.v3y = float(split[8])
        except: 
            self.v3y = 0.0
        
    def __getitem__(self, key):
        return self.__dict__[key]
    def write(self):
        string = "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n"%\
            (self.v1x,self.v2y,self.v3z,self.v1y,self.v1z,self.v2x,self.v2z
             ,self.v3x,self.v3y)
        return string
class Residue:
    """
    This will convert a list of Atom objects that belongs and will create a residue class taken into account the information of the first atom. It supoouse that all of them belong to the same residue of
    """
    def __init__(self, atoms):
        self.resName = atoms[0].resName
        self.resNumber = atoms[0].resNumber
        self.atoms = atoms

    def __getitem__(self, key):
        return self.__dict__[key]

class GRO:

    def __init__(self, file):
        self.file = file
        self.atoms = []
        self.name = ''
        self.NumberAtoms = 0
        self.vector = ""
        self.parse()

        
    def parse(self):
        with open(self.file, 'r') as f:
            lines = f.readlines()
        self.name = lines[0].strip()
        self.NumberAtoms = int(lines[1])
        self.vector = VECTOR(lines[-1])
        
        if len(lines) != self.NumberAtoms + 3: raise InputERROR('The number of atoms do not match or is something wrong with the gro file. See: https://manual.gromacs.org/archive/5.0.4/online/gro.html')
        
        for line in lines[2:-1]:
            atom = Atom(line)
            self.atoms.append(atom)
            
        
    def get_atoms(self, to_dict = True):
        """Return a list of all atoms.

        If to_dict is True, each atom is represented as a dictionary.
        Otherwise, a list of Atom objects is returned."""
        
        if to_dict: return [x.__dict__ for x in self.atoms]
        else: return self.atoms
        

    def get_residues(self):
        """Rearrange the data as residues function.

        Returns:
            list: a list of Residues objects
        """
        residues = []
        i = 0
        while i < len(self.atoms):
            resAtoms = []
            resName = self.atoms[i].resName
            resNumber = self.atoms[i].resNumber
            j = i
            while self.atoms[j].resName == resName and self.atoms[j].resNumber == resNumber:
                resAtoms.append(self.atoms[j])
                j += 1
                if j >= len(self.atoms): break
            residues.append(Residue(resAtoms))
            i = j
        
        return residues

    def edges(self):
        x_min = min(self.atoms, key=lambda item: item.x)
        y_min = min(self.atoms, key=lambda item: item.y)
        z_min = min(self.atoms, key=lambda item: item.z)

        x_max = max(self.atoms, key=lambda item: item.x)
        y_max = max(self.atoms, key=lambda item: item.y)
        z_max = max(self.atoms, key=lambda item: item.z)
        
        print(f"""
              X axes:
                     Length: {round(x_max.x-x_min.x, 3)} nm
                     Bounding atoms: {x_max.atomName}_{x_max.atomNumber} ({x_max.resName}) ----- {x_min.atomName}_{x_min.atomNumber} ({x_min.resName})
              Y axes:
                     Length: {round(y_max.y-y_min.y, 3)} nm
                     Bounding atoms: {y_max.atomName}_{y_max.atomNumber} ({y_max.resName}) ----- {y_min.atomName}_{y_min.atomNumber} ({y_min.resName})
              Z axes:
                     Length: {round(z_max.z-z_min.z, 3)} nm
                     Bounding atoms: {z_max.atomName}_{z_max.atomNumber} ({z_max.resName}) ----- {z_min.atomName}_{z_min.atomNumber} ({z_min.resName})
            """)
        
    def box_center(self, resName_selection = [], H_atoms = True):
        if resName_selection:
            atoms = [atom for atom in self.atoms if atom.resName in resName_selection]
        else:
            atoms = self.atoms
        
        if not H_atoms:
            atoms = [atom for atom in atoms if atom.atomName[0] not in 'Hh']

        x_min = min(atoms, key=lambda item: item.x)
        y_min = min(atoms, key=lambda item: item.y)
        z_min = min(atoms, key=lambda item: item.z)

        x_max = max(atoms, key=lambda item: item.x)
        y_max = max(atoms, key=lambda item: item.y)
        z_max = max(atoms, key=lambda item: item.z)
         
        return (round((x_max.x+x_min.x)/2,3), round((y_max.y+y_min.y)/2,3), round((z_max.z+z_min.z)/2,3))

    def centroid(self, resName_selection = [], H_atoms = True):
        
        if resName_selection:
            atoms = [atom for atom in self.atoms if atom.resName in resName_selection]
        else:
            atoms = self.atoms
        
        if not H_atoms:
            atoms = [atom for atom in atoms if atom.atomName[0] not in 'Hh']
        
        centroid = (sum([np.array([atom.x,atom.y,atom.z]) for atom in atoms]) / len(atoms))
        centroid = [round(coord, 3) for coord in centroid]
        return tuple(centroid)
    
    def write(self, out_file = None, backup = True):
        #This could be used if some modifications were made in the instance
        if out_file:
            to_work = out_file
        else:
            to_work = self.file
            if backup: tools.backoff(to_work)

        with open(to_work,'w') as out:
            if self.name: out.write(self.name+'\n')
            out.write(str(self.NumberAtoms)+'\n')
            for atom in self.atoms:
                out.write("%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n"%
                          (atom.resNumber,atom.resName,atom.atomName, 
                          atom.atomNumber,atom.x,atom.y,atom.z,atom.vX,atom.vY,
                          atom.vZ
                          ))            
            if self.vector: out.write(self.vector.write())
            
            
if __name__ == '__main__':
    #import sys
    #gro = GRO(sys.argv[1])
    gro = GRO('confout.gro')
    atoms = gro.get_atoms(to_dict = False)
    gro.edges()
    print(gro.box_center())
    
    for atom in gro.atoms:
        if atom.resName == 'LIA':
            atom.x = 0
            atom.y = 0
            atom.z = 0

    
    gro.write()
    #print(gro.NumberAtoms)
    #print(gro.name)
   
    



 
