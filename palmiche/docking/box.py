#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import palmiche.utils.pdb as pdb
from palmiche.utils.assembleMD import fix_pdb
import palmiche.utils.tools as tools
import os.path




def databox(input_file, exp_structure, pymol_path_cmd = "pymol", lig_name = None, distance = 0.4):

    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"The file {input_file} doesn't exist or is not accesible")
    if not os.path.isfile(exp_structure):
        raise FileNotFoundError(f"The file {exp_structure} doesn't exist or is not accesible")

    """
    This function (fix_pdb) will export the pdb file:
        1-PROTEIN.pdb
    The next steps on this function will make use of PROTEIN.pdb

    """
    PDB = fix_pdb(input_file)
    tools.run(f"{pymol_path_cmd} PROTEIN.pdb {exp_structure} &")
    answer1 = "a"
    while answer1 != "c":
        answer1 = input("""You need to align the pdb structure of the PDB file that
          has the ligands to the structure for fix_pdb(input_file),
          and create and save a
          pdb file under the name PROTEIN.pdb (the same name as the
                exported structure) that has the protein and the ligands
          aligned. This could be done in Pymol.
          0-Open input_file and the experimental pdb. In this order, so, the ligands
           will be at the end of the file. (this should had done
            automathically)
          1-Align the experimental or the one that has the ligands orientation
          to input_file, center the view in PROTEIN to seee the alignment.
          2-Slect all the secuence of input_file
          3-Select also the ligands
          4-Export the selection as input_file (Replacing the existing PDB file)\n\n

          If you already finished this step type 'c' to continue: """).lower()

    """
    Rectified again the PDB is needed becasue Pymol change the position of the
    TER flag.
    """

    PDB = pdb.PDB("PROTEIN.pdb")#I nead to load the new PDB from the pymol work
    PDB.fix_chains()

    #Now the name of the ligands will be changed
    if lig_name:
        for atom in PDB.atoms:
            if atom.resName == lig_name:
                atom.resName = f"LI{atom.chainID}"
    #write the new PDB file with the Ligan names corrected
    PDB.write()
    chains = PDB.get_chains()
    #Creating the index file and  the pdb to get the box dimensions
    #The first step is the generations of the chain index with make_ndx
    with open("opt", "w") as opt:
        for chain in chains:
            opt.write(f"chain {chain}\n")
        opt.write("q\n")
    tools.run("cat opt | gmx make_ndx -f PROTEIN.pdb -o index")

    with open("opt", "w") as opt:
        for chain in chains:
            opt.write(f"\"SELECT_{chain}\" group \"ch{chain}\" and same residue as within {distance} of resname LI{chain};\n")
    tools.run("gmx select -s PROTEIN.pdb -sf opt -n index.ndx -on index.ndx")

    #deleting the line _f0_t0.000 in the file
    with open("index.ndx", "r") as index:
        data = index.read()
        data = data.replace("_f0_t0.000","")
    with open("index.ndx", "w") as index:
        index.write(data)


    #create the PDB selection, get the center of the box and edges, delete the pdb file
    boxes = {}
    for chain in chains:
        tools.run(f"echo \"SELECT_{chain}\" | gmx editconf -f PROTEIN.pdb -o SELECT_{chain}.pdb -n index.ndx")
        SELECT = pdb.PDB(f"SELECT_{chain}.pdb")
        boxes[chain] = {"boxcenter": SELECT.box_center(), "boxsize": SELECT.edges()}
        tools.rm(f"SELECT_{chain}.pdb")
    [tools.rm(item) for item in ["opt", "index*", "#*"]]

    return boxes

if __name__ == '__main__':
    boxes = databox("step5_input.gro", "7e27.pdb", pymol_path_cmd = "/home/ale/pymol/pymol", lig_name = "HV6")
    print(boxes)
