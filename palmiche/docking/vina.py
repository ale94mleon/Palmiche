#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob as glob
import multiprocessing as mp
import os
import numpy as np
from palmiche.utils import tools
import tqdm

#========Deffinitions of class and methosd to get the output of Vina===========
class Atom:
    #https://userguide.mdanalysis.org/stable/formats/reference/pdbqt.html#writing-out
    def __init__(self, line):
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

class CHUNK_VINA_OUT:
    def __init__(self, chunk):
        self.chunk = chunk
        self.atoms = []
        self.run = None
        self.freeEnergy = None
        self.RMSD1 = None
        self.RMSD2 = None
        self.parse()

    def parse(self):
        for line in self.chunk:
            if line.startswith("MODEL"):
                self.run = int(line[5:])
            elif line.startswith("REMARK VINA RESULT:"):
                    (self.freeEnergy, self.RMSD1, self.RMSD2) = [float(number) for number in line.split(":")[-1].split()]
                    
            elif line.startswith("ATOM"):
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
                f.writelines(self.chunk)
        else:
            with open(f"Run_{self.run}.pdbqt","w") as f:
                f.writelines(self.chunk)            

class VINA_OUT:
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

            if lines[i].startswith("MODEL"):
                j = i
                tmp_chunk = []
                while (not lines[j].startswith("ENDMDL")) and (j < len(lines)):
                    tmp_chunk.append(lines[j])
                    j += 1
                    i += 1
                tmp_chunk.append("ENDMDL\n")
                self.chunks.append(CHUNK_VINA_OUT(tmp_chunk))

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
        PO_chunk and NO_chunk obkjects and positive_oriented.pdbqt and negative_oriented.pdbq

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
            PO_chunk = min(positive_orientation,key=lambda x: x.freeEnergy)
            PO_chunk.write("positive_oriented.pdbqt")
        else:
            print("A positive orientaion was not found")
        if negative_orientation:    
            NO_chunk = min(negative_orientation,key=lambda x: x.freeEnergy)
            NO_chunk.write("negative_oriented.pdbqt")
        else:
            print("A negative orientaion was not found")
        
        return PO_chunk, NO_chunk
            
    def BestEnergy(self, write = True):
        min_chunk = min(self.chunks, key= lambda x: x.freeEnergy)
        if write: min_chunk.write("best_energy.pdbqt")
        
        return min_chunk
            
        
#==============================================================================





#====================Congfiguration============================================
def vinaconfig(receptor_path, ligand_path, grid_center, size, modes = 1, config_path = "./", chain = "A"):
    """
    

    Parameters
    ----------
     receptor : TYPE
        DESCRIPTION:
            path to the receptor (relative or absolute)
    
    ligand : TYPE
        DESCRIPTION:
            path to the ligand (relative or absolute)
            e.g:
                mol.pdbqts

    grid_center : TYPE
        DESCRIPTION:
            A tupla with the float coords of the center:
                (x,y,z)
    size : TYPE
        DESCRIPTION:
            A tupla with the size (float > 0) on each coord:
                (Sx,Sy,Sz)
            It is important to consider that the grid space used is grid_space = 0.375
            This is the default value used 
    modes : TYPE, optional
        DESCRIPTION:
            How many modes will be computed
            The default is 1. 

    Returns
    -------
    None.

    """


    name = os.path.basename(os.path.abspath(ligand_path)).split(".")[0]
    config_path = os.path.abspath(config_path)
    receptor_path = os.path.abspath(receptor_path)
    ligand_path = os.path.abspath(ligand_path)
    
    final=open(f"{config_path}/{name}_{chain}.config",'w')

    final.write(f"receptor      =  {receptor_path}\n")
    final.write(f"ligand        =  {ligand_path}\n")
    
    final.write(f"center_x      =  {grid_center[0]}\n")
    final.write(f"center_y      =  {grid_center[1]}\n")
    final.write(f"center_z      =  {grid_center[2]}\n")
    
    final.write(f"size_x        =  {size[0]}\n")
    final.write(f"size_y        =  {size[1]}\n")
    final.write(f"size_z        =  {size[2]}\n")
    final.write(f"num_modes     =  {modes}\n")
    
    final.write(f"out       =  {config_path}/{name}_{chain}_out.pdbqt\n")

    final.close()


#==============================================================================

#=================Deffinition of the run=======================================
def runvina(grid_opt, exhaustiveness = 25, cpu = 4, chains = ["A"], modes = 1, continue_check = False):
    """
    Parameters
    ----------
    grid_opt : TYPE
        DESCRIPTION.
    exhaustiveness : TYPE, optional
        DESCRIPTION. The default is 25.
    cpu : TYPE, optional
        DESCRIPTION. The default is 4.
    chains : TYPE, optional
        DESCRIPTION. The default is ["A"].

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    cwd = os.getcwd()
    
    
    ligand_paths = glob.glob(os.path.join(cwd, "ligands", "*.pdbqt"))
    ligands_name = [os.path.basename(ligand_path).split(".")[0] for ligand_path in ligand_paths]
    
    receptor_paths = glob.glob(os.path.join(cwd, "receptors", "*.pdbqt"))
    receptors_name = [os.path.basename(receptor_path).split(".")[0] for receptor_path in receptor_paths]
    
    cmd = []
    for i, receptor_name in enumerate(receptors_name):
        
        tools.makedirs(f"{cwd}/run_vina_run/{receptor_name}")
        for j, ligand_name in enumerate(ligands_name):
            tools.makedirs(f"{cwd}/run_vina_run/{receptor_name}/{ligand_name}")
            for chain in chains:
                tools.makedirs(f"{cwd}/run_vina_run/{receptor_name}/{ligand_name}/{chain}")
                #print(f"{cwd}/run_vina_run/{receptor_name}/{ligand_name}/{chain}")
                #tools.cp(receptor_paths[i], f"{cwd}/run_vina_run/{receptor_name}/{ligand_name}/{chain}")
                #tools.cp(ligand_paths[j], f"{cwd}/run_vina_run/{receptor_name}/{ligand_name}/{chain}")
                try:
                    param = grid_opt[receptor_name][chain]
                    
                    vinaconfig(receptor_paths[i],ligand_paths[j],
                                               param["boxcenter"],
                                               param["boxsize"],
                                               config_path=f"{cwd}/run_vina_run/{receptor_name}/{ligand_name}/{chain}",
                                               chain=chain,
                                               modes=modes)
                    #Check for the _out.pdbqt, this is only exported when the calculation finished
                    if continue_check:
                        if not os.path.isfile(f"{cwd}/run_vina_run/{receptor_name}/{ligand_name}/{chain}/{ligand_name}_{chain}_out.pdbqt"):
                            cmd.append(f"vina --config {cwd}/run_vina_run/{receptor_name}/{ligand_name}/{chain}/{ligand_name}_{chain}.config --exhaustiveness={exhaustiveness} --cpu={cpu}")
                    else:
                        cmd.append(f"vina --config {cwd}/run_vina_run/{receptor_name}/{ligand_name}/{chain}/{ligand_name}_{chain}.config --exhaustiveness={exhaustiveness} --cpu={cpu}")
                except:
                    print(f"{receptor_name} or/and {chain} of {receptor_name} is/are not present(s) in {grid_opt}")
    
    
    #Parallel running
    #Creating the iterable was already created (cmd)
    #[print(c) for c in cmd]
    #Creating the pool
    pool = mp.Pool(int(mp.cpu_count()/cpu))
    #Running
    for i in tqdm.tqdm(pool.imap_unordered(tools.run, cmd), total=len(cmd)):
        pass
    
    #Closing
    pool.close()
#==============================================================================



if __name__ == '__main__':
    pass
    #vina = VINA_OUT("mol1_out.pdbqt")
    #vina.PosNegConf(atom_1 = 17 , atom_2 = 1)
    #runvina("/home/ale/MY_PYTHON_PACKEGES/palmiche/examples/Vina_Docking/ref_info/boxes.box")
