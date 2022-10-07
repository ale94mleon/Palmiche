#!/usr/bin/env python

from palmiche.utils import tools
from string import ascii_uppercase
import numpy as np
"""
Taken from https://charmm-gui.org/?doc=lecture&module=scientific&lesson=7
"""

class Atom:
    """
    https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM    
    """
    def __init__(self, line):
        self.code = "ATOM"
        self.serial = int(line[6:11])
        self.name = line[12:16].strip()
        self.altLoc = line[16]
        self.resName = line[17:21].strip()
        self.chainID = line[21]
        self.resSeq = int(line[22:26])
        self.iCode = line[26]
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occupancy = line[54:60].strip()
        self.tempFactor = line[60:66].strip()
        self.element = line[76:78].strip()
        self.charge = line[78:80].strip()
        if self.occupancy: self.occupancy = float(self.occupancy)
        if self.tempFactor: self.tempFactor = float(self.tempFactor)

    def __getitem__(self, key):
        return self.__dict__[key]


class Hetatm:
    """
    https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM    
    """
    def __init__(self, line):
        self.code = "HETATM"
        self.serial = int(line[6:11])
        self.name = line[12:16].strip()
        self.altLoc = line[16]
        self.resName = line[17:21].strip()
        self.chainID = line[21]
        self.resSeq = int(line[22:26])
        self.iCode = line[26]
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occupancy = line[54:60].strip()
        self.tempFactor = line[60:66].strip()
        self.element = line[76:78].strip()
        self.charge = line[78:80].strip()
        if self.occupancy: self.occupancy = float(self.occupancy)
        if self.tempFactor: self.tempFactor = float(self.tempFactor)

    def __getitem__(self, key):
        return self.__dict__[key]


class CRYST1:
    """
    https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html#CRYST1    
    """
    def __init__(self, line):
        self.a = float(line[6:15])			    #Real(9.3)     a              a (Angstroms).
        self.b = float(line[15:24])			    #Real(9.3)     b              b (Angstroms).
        self.c = float(line[24:33])			    #Real(9.3)     c              c (Angstroms).
        self.alpha = float(line[33:40])			#Real(7.2)     alpha          alpha (degrees).
        self.beta = float(line[40:47])			#Real(7.2)     beta           beta (degrees).
        self.gamma = float(line[47:54])			#Real(7.2)     gamma          gamma (degrees).
        self.sGroup = line[55:66]			    #LString       sGroup         Space  group.
        try:
            self.z = int(line[66:70])			    #Integer       z              Z value.
        except:
               self.z = "" 

    def __getitem__(self, key):
        return self.__dict__[key]
    def write(self):
        string = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%-12s%4s\n"%\
            (self.a,self.b,self.c,self.alpha,self.beta,self.gamma,self.sGroup,self.z)
        return string
class Residue:
    """
    This will convert a list of Atom objects that belongs and will create a residue class taken into account the information of the first atom. It supoouse that all of them belong to the same residue of
    """
    def __init__(self, atoms):
        self.resName = atoms[0].resName
        self.resSeq = atoms[0].resSeq
        self.chainID = atoms[0].chainID
        self.atoms = atoms

    def __getitem__(self, key):
        return self.__dict__[key]

class PDB:
    def __init__(self, file):
        self.file = file
        self.title = ""
        self.cryst1 = CRYST1("0"*54+"  P 1 21 1     1")#Just to initialized to a generci cell, see the web referenced in CRYST1
        self.atoms = []
        self.parse()


    def parse(self):
        MODEL = None
        
        f = open(self.file, 'r')
        
        for line in f.readlines():

            if line.startswith('TITLE'): self.title = line.strip()
            if line.startswith('CRYST1'): self.cryst1 = CRYST1(line)
            if line.startswith('MODEL'): MODEL = int(line.split()[1])
            
            if line.startswith('ATOM'):
                atom = Atom(line)
                atom.MODEL = MODEL
                self.atoms.append(atom)
            if line.startswith('HETATM'):
                atom = Hetatm(line)
                atom.MODEL = MODEL
                self.atoms.append(atom)
        f.close()
        


    def get_chains(self):
        chains = set()
        for atom in self.atoms:
            chains.add(atom.chainID)
        return sorted(chains)#This atomathically convert to list
        

    def get_atoms(self):
        """Return a list of all atoms.
        each atom is represented as a dictionary.
        """
        return [x.__dict__ for x in self.atoms]
        

    def get_model(self, model_num, to_dict=True):
        """Return all atoms where MODEL == model_num"""
        model_atoms = [x for x in self.atoms if x.MODEL == model_num]
        if to_dict:
            return [atom.__dict__ for atom in model_atoms]
        else:
            return model_atoms
    

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
            resSeq = self.atoms[i].resSeq
            chainID = self.atoms[i].chainID
            j = i
            while self.atoms[j].resName == resName and self.atoms[j].resSeq == resSeq and self.atoms[j].chainID == chainID:
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
        
        length = [round(x_max.x-x_min.x, 3), round(y_max.y-y_min.y, 3), round(z_max.z-z_min.z, 3)]
        
        print(f"""
X axes:
       Length: {length[0]} Angstroms
       Bounding atoms: {x_max.name}_{x_max.serial} ({x_max.resName}_{x_max.resSeq}-{x_max.chainID}) ----- {x_min.name}_{x_min.serial} ({x_min.resName}_{x_min.resSeq}-{x_min.chainID})
Y axes:
       Length: {length[1]} Angstroms
       Bounding atoms: {y_max.name}_{y_max.serial} ({y_max.resName}_{y_max.resSeq}-{y_max.chainID}) ----- {y_min.name}_{y_min.serial} ({y_min.resName}_{y_min.resSeq}-{y_min.chainID})
Z axes:
       Length: {length[2]} Angstroms
       Bounding atoms: {z_max.name}_{z_max.serial} ({z_max.resName}_{z_max.resSeq}-{y_max.chainID}) ----- {z_min.name}_{z_min.serial} ({z_min.resName}_{z_min.resSeq}-{z_min.chainID})
            """)
        return tuple(length)
        
    def box_center(self, resName_selection = [], H_atoms = True):
        if resName_selection:
            atoms = [atom for atom in self.atoms if atom.resName in resName_selection]
        else:
            atoms = self.atoms
        
        if not H_atoms:
            atoms = [atom for atom in atoms if atom.name[0] not in 'Hh']
        
        x_min = min(atoms, key=lambda item: item.x)
        y_min = min(atoms, key=lambda item: item.y)
        z_min = min(atoms, key=lambda item: item.z)

        x_max = max(atoms, key=lambda item: item.x)
        y_max = max(atoms, key=lambda item: item.y)
        z_max = max(atoms, key=lambda item: item.z)
         
        return round((x_max.x+x_min.x)/2,3), round((y_max.y+y_min.y)/2,3), round((z_max.z+z_min.z)/2,3)
    
    def centroid(self, resName_selection = [], H_atoms = True):
        if resName_selection:
            atoms = [atom for atom in self.atoms if atom.resName in resName_selection]
        else:
            atoms = self.atoms
        
        if not H_atoms:
            atoms = [atom for atom in atoms if atom.name[0] not in 'Hh']
        
        centroid = (sum([np.array([atom.x,atom.y,atom.z]) for atom in atoms]) / len(atoms))
        centroid = [round(coord, 3) for coord in centroid]
        return tuple(centroid)


    def fix_chains(self, HETATM = False):
        #In the case that you ahve only one aa per chain with the same residue
        #number this will fail. It will asign chain A to all. 
        
        #Here we are using only the ATOMS for the CHAIN rectificatin, not the HETATM
        #If HETATM False, use only ATOMS,
        #If HETATM True, use ATOMS and HEATM
        if not HETATM:
            atoms = [atom for atom in self.atoms if atom.code == "ATOM"]
        else:
            atoms = self.atoms
        
        chain = 0    
        for (i, atom) in enumerate(atoms):
            if i == 0:
                atom.chainID = ascii_uppercase[chain]
            elif atoms[i].resSeq >= atoms[i-1].resSeq:
                atom.chainID = ascii_uppercase[chain]
            else:
                chain += 1
            atom.chainID = ascii_uppercase[chain]
    
    def add_element(self):
        """
        This function add the columnd elemnt to the pdb data

        Returns
        -------
        None.

        """
        for atom in self.atoms:
            if not atom.element:
                element = ""
                for char in atom.name:
                    if char.isalpha():
                        element = char.upper()
                        break
                atom.element = element
    

    def write(self, out_file = None, backup = True):
        #This could be used if some modifications were made in the instance
        if out_file:
            to_work = out_file
        else:
            to_work = self.file
            if backup: tools.backoff(to_work)
        with open(to_work,'w') as out:
            if self.title: out.write(self.title+'\n')
            if self.cryst1.a or self.cryst1.b or self.cryst1.c: #If the length of at least one vector is greater than 0
                out.write(self.cryst1.write())
            if self.atoms[0].MODEL != None: out.write("MODEL %6i\n"%(self.atoms[0].MODEL))
            
            
            """
            For the case that the protein have more than one chain, teh TER 
            flag will be printed. If not (thae protein has only one chain or 
            the chain is not decleared) the TER flag will not printed
            """
            CHAINS_ID = self.get_chains()#Becasue get_chains belongs to the class you need to call it like a self.function(), the method of the class
            if len(CHAINS_ID) > 1:
                last = self.atoms[0].chainID
                for (i,atom) in enumerate(self.atoms):
                    
                    #INCREDIBLE, PropKa if the atoms name is shifted to the left dont understand!!!
                    #For that reason I need to change when the len is 3. NO IDEA!!!!
                    if len(atom.name) == 3:
                        string_name = f"{atom.name:>4}"
                    else:
                        string_name = f"{atom.name:^4}"
                        
                    if atom.chainID != last:
                        out.write("TER\n")
                        last = atom.chainID

                    #This check is necesary for the case that some residues has a name greater than 3 letters (ASPH)
                    #This generate something like ASPHA, where A is the Chain identifier
                    
                    out.write(f"{atom.code:<6}{atom.serial:5} {string_name}"
                              f"{atom.altLoc:1}{atom.resName:<4}{atom.chainID:1}"
                              f"{atom.resSeq:4}{atom.iCode:1}   {atom.x:8.3f}"
                              f"{atom.y:8.3f}{atom.z:8.3f}{atom.occupancy:6.2f}"
                              f"{atom.tempFactor:6.2f}           {atom.element:2}"
                              f"{atom.charge:2}\n")
                out.write("TER\n")#This is for the last TER flag
            else:
                for (i,atom) in enumerate(self.atoms):
                    #INCREDIBLE, PropKa if the atoms name is shifted to the left dont understand!!!
                    if len(atom.name) == 3:
                        string_name = f"{atom.name:>4}"
                    else:
                        string_name = f"{atom.name:^4}"
                        
                    out.write(f"{atom.code:<6}{atom.serial:5} {string_name}"
                              f"{atom.altLoc:1}{atom.resName:<4}{atom.chainID:1}"
                              f"{atom.resSeq:4}{atom.iCode:1}   {atom.x:8.3f}"
                              f"{atom.y:8.3f}{atom.z:8.3f}{atom.occupancy:6.2f}"
                              f"{atom.tempFactor:6.2f}           {atom.element:2}"
                              f"{atom.charge:2}\n")
                        

if __name__ == '__main__':
    #import sys
    #pdb = PDB(sys.argv[1])
    #pdb = PDB('PROTEIN.pdb')
    #print(pdb)
    #atoms = pdb.atoms
    #print(atoms[0].__dict__)
    #print(pdb.get_chains())
    #for atom in atoms: #Could be used to rectificate the pdb and put the chain information
    #    atom.chainID = "G" 
    #print(pdb.get_chains())
    #pdb.cryst1.a = 500
    #print(pdb.cryst1.write())
    #for atom in atoms:
    #    if atom.resName == "LIA":
    #        atom.x = 0
    #        atom.y = 0
    #        atom.z = 0
    #pdb.write("j.pdb")
    #print(isinstance(pdb, PDB))
    #print(type(pdb).__name__)
    

    #print the name of each atom
    #for atom in atoms:
    #    print(atom.serial, atom.name, atom.x, atom.y, atom.z, atom.element)
    #lig_name = "HV6"
    #for atom in pdb.atoms:
    #    if atom.resName == lig_name:
    #        atom.resName = f"LI{atom.chainID}"
    #print(pdb.centroid(), pdb.box_center(), pdb.edges())
    #rint(pdb.atoms[-1].__dict__)
    pass
