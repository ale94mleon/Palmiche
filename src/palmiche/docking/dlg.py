#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
class Atom:
    #https://userguide.mdanalysis.org/stable/formats/reference/pdbqt.html#writing-out
    def __init__(self, line):
        if "DOCKED: " in line: line = line[len("DOCKED: "):]
        self.lineType = "ATOM"
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
        self.partialChrg = line[66:76].strip()
        self.atomType = line[78:80].strip()
    def __getitem__(self, key):
        return self.__dict__[key]

class Hetatm:
    def __init__(self, line):
        if "DOCKED: " in line: line = line[len("DOCKED: "):]
        self.lineType = "HETATM"
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
        self.partialChrg = line[66:76].strip()
        self.atomType = line[78:80].strip()
    def __getitem__(self, key):
        return self.__dict__[key]

class Remark:
    def __init__(self, line):
        pass

class CHUNK_DLG:
    def __init__(self, chunk):
        self.chunk = chunk
        self.atoms = []
        self.run = None
        self.freeEnergy = None
        self.to_write = ""
        self.parse()

    def parse(self):
        slip = len("DOCKED: ")
        for line in self.chunk:
            self.to_write += line[slip:]

            if line.startswith("DOCKED: USER"):
                if line[16:19] == "Run":
                    self.run = int(line[21:])
                elif "Estimated Free Energy of Binding" in line:
                    self.freeEnergy = float(line[53:61])
                else:
                    pass
            elif line.startswith("DOCKED: ATOM"):
                self.atoms.append(Atom(line))
            else:
                pass

    def get_atoms(self):
        """Return a list of all atoms.

        If to_dict is True, each atom is represented as a dictionary.
        Otherwise, a list of Atom objects is returned."""
        return [x.__dict__ for x in self.atoms]

    def get_vector(self, atom_index1, atom_index2, normalized = True):
        """
        Parameters
        ----------
        atom_index1 : TYPE: int
            DESCRIPTION: The index of the first atom. Remember that python list
            start at 0. So, the atom number 65 in a pdbqt will be at the position
            64 of the list of atoms.
        atom_index2 : TYPE: int
            DESCRIPTION: The index of the second atom. Remember that python list
            start at 0. So, the atom number 1 in a pdbqt will be at the position
            0 of the list of atoms.
        normalized : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        a numpy array. the vector between the two atoms, the vector point to 
        to atom2 (atom1|----->atom2)
        """
        XYZ_atom1 = (self.atoms[atom_index1].x, self.atoms[atom_index1].y, self.atoms[atom_index1].z)
        XYZ_atom2 = (self.atoms[atom_index2].x, self.atoms[atom_index2].y, self.atoms[atom_index2].z)
        vector = np.array(XYZ_atom2) - np.array(XYZ_atom1)
        
        if normalized:
            vector /= np.linalg.norm(vector)
        return vector
        
    def write(self, name = None):
        if name:
            with open(name,"w") as f:
                f.write(self.to_write)
        else:
            with open(f"Run_{self.run}.pdbqt","w") as f:
                f.write(self.to_write)            

class DLG():
    """
    To acces the chunks you need to take into account that 
    """
    def __init__(self, file):
        self.file = file

        self.chunks = []
        self.parse()

    def parse(self):
        with open(self.file, "r") as input_file:
            lines = input_file.readlines()
        i = 0
        while i < len(lines):

            if "FINAL GENETIC ALGORITHM DOCKED STATE" in lines[i]:
                j = i
                tmp_chunk = []
                while (not lines[j].startswith("_")) and (j < len(lines)):

                    if lines[j].startswith("DOCKED:"):
                        tmp_chunk.append(lines[j])
                    j += 1
                    i += 1
                self.chunks.append(CHUNK_DLG(tmp_chunk))

            i += 1
    def PosNegConf(self, atom_1, atom_2, vector = (0,0,1)):
        """
        

        Parameters
        ----------
        atom_1 : TYPE: int
            DESCRIPTION: The index of the atom 1 (you could know this using pymol, label>atomidentifier>ID )
        atom_2 : TYPE: int
            DESCRIPTION: The index of the atom 2 (you could know this using pymol, label>atomidentifier>ID )
        vector : TYPE, optional, is the reference vector in order to 
        look for the confomers in the direction and in the opposite direction 
        of this vector
            DESCRIPTION. The default is [0,0,1].

        Returns
        -------
        positive_oriented.pdbqt
        negative_oriented.pdbq

        """
        
        positive_orientation = []
        negative_orientation = []
     
        
        for chunk in self.chunks:
            vector_chunk = chunk.get_vector(atom_1 - 1, atom_2 - 1, normalized=True)
            dot_product = np.dot(vector_chunk, vector)
            
            if dot_product >= 0:
                positive_orientation.append(chunk)
            else:
                negative_orientation.append(chunk)

        if positive_orientation:
            min(positive_orientation,key=lambda x: x.freeEnergy).write("positive_oriented.pdbqt")
        else:
            print("A positive orientaion was not found")
        if negative_orientation:    
            min(negative_orientation,key=lambda x: x.freeEnergy).write("negative_oriented.pdbqt")
        else:
            print("A negative orientaion was not found")          

if __name__ == '__main__':
    oriented_vector = np.array([0,0,1])
    dlg = DLG("cis1_bh267_Cell_682.dlg")
    dlg.PosNegConf(atom_1 = 17 , atom_2 = 1)
                    