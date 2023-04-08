#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from palmiche.utils import tools
from palmiche.utils import pdb
from palmiche.utils.pKa2pdb2gmx import pka2gmx
from glob import glob
import tempfile



def get_path_dict(ligands_path = "./ligands", receptors_path = "./receptors", dockings_path = "./run_vina_run"):
    
    ligands_path = os.path.abspath(ligands_path)
    receptors_path = os.path.abspath(receptors_path)
    dockings_path = os.path.abspath(dockings_path)
    
    """These assume that:
        in receptor path the tree is:
            receptor1/
            receptor2/
            receptor3/
            (NO more directories)
        And inside of this directories there are three files:
            receptor1.pdb, step5_input.gro, step5_input.pdb and MEMBRANE.pdb
            
        As you see, this assumes that the receptor1.pdb has the same name as
        the parental directory (receptor1).
        
        in ligand path the tree is:
            ligand1/
            ligand2/
            ligand3/
            (NO more directories)
        And inside of this directories there are four files
        list.lis, is the list to change the order of the ligand, this is only the path, 
            In further step this file need to be created
            mol.frcmod
            mol.mol2
            mol.pdb
            list.lis
        
        in docking path the tree is the one generated with the script in palmiche/examples/Vina_Docking
        receptor1
            ligandX
                chainA
                    result1
                    .
                .   
            .    
        .
        """
    receptor_dict = {}
    for directory in tools.list_if_dir(receptors_path):
        
        receptor_dict[directory] = {"protein": os.path.join(receptors_path, directory, f"{directory}.pdb"),
                                    "step5_input_gro": os.path.join(receptors_path, directory, "step5_input.gro"),
                                    "step5_input_pdb": os.path.join(receptors_path, directory, "step5_input.pdb"),
                                    "membrane": os.path.join(receptors_path, directory, "MEMBRANE.pdb")}
    
   
    ligand_dict = {}
    for directory in tools.list_if_dir(ligands_path):
        ligand_dict[directory] = {"frcmod": os.path.join(ligands_path, directory, "mol.frcmod"),
                                  "mol2": os.path.join(ligands_path, directory, "mol.mol2"),
                                  "pdb": os.path.join(ligands_path, directory, "mol.pdb"),
                                  "list": os.path.join(ligands_path, directory, "list.lis")}
    
    docking_dict = {}
    for receptor in tools.list_if_dir(dockings_path):
        for ligand in tools.list_if_dir(os.path.join(dockings_path, receptor)): #In this way if the receptor where teasted to different ligands  (een so, they must be in the deffinition of the ligands_dict) this will tak it intop account
            for chain in tools.list_if_dir(os.path.join(dockings_path, receptor, ligand)):
                #Creating nested dictionary
                docking_dict.setdefault(receptor,{}).setdefault(ligand,{})[chain] = {"BE":os.path.join(dockings_path, receptor, ligand, chain, f"{ligand}_{chain}_BE.pdbqt"),
                                                                                     "NO":os.path.join(dockings_path, receptor, ligand, chain, f"{ligand}_{chain}_NO.pdbqt"),
                                                                                     "PO":os.path.join(dockings_path, receptor, ligand, chain, f"{ligand}_{chain}_PO.pdbqt")
                                                                                     }
    return receptor_dict, ligand_dict, docking_dict

def fix_pdb(step5_input_pdb, out = "PROTEIN.pdb"):
    """
    Take a pdb file and return the protein with the chains rectified, also convert CHARMM's atoms name to AMBER
    Parameters
    ----------
    step5_input_gro :  name of the file or path to the configuration file
        DESCRIPTION.A configuration file exported by CHARMM-GUI, gro or pdb.
        It is useful to use the last gro file because this one has the information
        for the construction of the gmx box on further steps.
    out : The name or path toe the output rectified pdb
 
    Returns
    -------
    1- PDB_OBJECT, this is a PDB object
    2- Also, a file called PROTEIN.pdb will be exported,
    This pdb file is protonated, with the chain rectified and also the element column added
    """
    if not os.path.isfile(step5_input_pdb):
        raise FileNotFoundError(f"The file {step5_input_pdb} doesn't exist or is not accessible")
    
    tmppdb1 = tempfile.NamedTemporaryFile(suffix='.pdb')
    tmppdb2 = tempfile.NamedTemporaryFile(suffix='.pdb')
    tmpndx = tempfile.NamedTemporaryFile(suffix='.ndx')

    ext = step5_input_pdb.split(".")[-1]
    if ext == "gro":
        tools.run(f"export GMX_MAXBACKUP=-1; echo gmx editconf -f {step5_input_pdb} -o {tmppdb1.name}")
        
        with open("warning.txt", "w") as w:
            w.write(f"""
            WARNING! 
            This is a warning from the function fix_pdb of assembleMD module
            The file {step5_input_pdb} is a gro file and
            gmx editconf was used to convert to pdb. If this function was used
            to start the assemble from CHARMM-GUI, please select the
            step5_input.pdb. Because the gro file doesn't have the information
            related with the chains and is not possible to fix because the residues
            number is changed in CHARMM-GUI to a continues list and doesn't use
            the real residues number of the pdb"
            """)
    else:
        tools.cp(step5_input_pdb, tmppdb1.name)
    
    #Keeping only the protein
    tools.run(f"""
            export GMX_MAXBACKUP=-1
            echo "q" | gmx make_ndx -f {tmppdb1.name} -o {tmpndx.name}
            echo "Protein" | gmx editconf -f {tmppdb1.name} -o {tmppdb2.name} -n {tmpndx.name}
          """)
    
    #************Rectify atom names, chains, residue...************************
    #This are the translations dict the kkey are the charmm nomenclature and the
    #values are the amber nomenclature for the atoms of the patching groups
    #ACE: acetylated N-terminus
    #CT3_NME: methylamidated C-terminus
    #The lst dict is for rename the HIS resnames, the former will change the atoms,
    #the residue and also the residue number
    #I also need to repair the chain
    ACE_ACE = {"CAY":"CH3", "HY1":"HH31", "HY2":"HH32", "HY3":"HH33", "CY":"C", "OY":"O"}
    CT3_NME = {"NT":"N", "HNT":"H", "CAT":"CH3", "HT1":"HH31", "HT2":"HH32", "HT3":"HH33"}
    HSD_HIS = {"HSD":"HIS"}
    
    
    PDB_OBJECT = pdb.PDB(tmppdb2.name)
    PDB_OBJECT.fix_chains()
    PDB_OBJECT.add_element()
    
    atoms = PDB_OBJECT.atoms
    #Get the Chains
    
    for atom in atoms:
        
    
        if atom.name in ACE_ACE:
            atom.name = ACE_ACE[atom.name]
            atom.resName = "ACE"
            #-1 Because this is for the N terminus that start in the PDB, this could
            #change, take a look to your system. This is not completely optimal. 
            #In the case that the aa start at 0you will get a -1 of residue number. 
            #However, why to cap if you have th first residue?
            #There is not need to change the name of the residue
            atom.resSeq -= 1
        elif atom.name in CT3_NME:
            atom.name = CT3_NME[atom.name]
            atom.resName = "NME"
            atom.resSeq += 1
        elif atom.resName in HSD_HIS:
            atom.resName = HSD_HIS[atom.resName]
        else:
            pass
        
    #This automatically convert to list
    PDB_OBJECT.write(tmppdb2.name, backup = False)
    #Get correct protonation, here I use amber99sb-ildn because it is by default in GROMACS
    
    #Probably a better idea is use CHARMM force fields in irder to dont worry about the atom names, but this use (amber) is consistent with my project 
    pka2gmx(tmppdb2.name,
            output = tmppdb2.name, 
            protein_forcefield = "amber99sb-ildn",
            water_forcefield = "tip3p",
            ph = 7.0,
            pKa = "propka3",
            pKa_file_out = False,
            itp_top_out = False,
            num_file = False)
    PDB_OBJECT = pdb.PDB(tmppdb2.name) #Need to load the new pdb

    PDB_OBJECT.write(out, backup = False)
    
    return PDB_OBJECT

def get_cryst1(step5_input_gro_path):
    tools.run(f"""
              gmx editconf -f {step5_input_gro_path} -o tmp.pdb
              """)
    TMP_PDB = pdb.PDB("tmp.pdb")
    #I need to dived by 10 because the pdb is in Angstrom
    box_vector = (TMP_PDB.cryst1.a/10, TMP_PDB.cryst1.b/10, TMP_PDB.cryst1.c/10)
    box_angles = (TMP_PDB.cryst1.alpha, TMP_PDB.cryst1.beta, TMP_PDB.cryst1.gamma)
    tools.rm("tmp.pdb")
    return box_vector, box_angles
    
    
def gmxmemb(step5_input_pdb_path,
            membrane_out_path,
            membrane_select = "(group System and ! group TIP3) and ! group Protein",
            lipid_ff_path = '/home/ale/MY_PYTHON_PACKEGES/palmiche/data/GROMACS.ff/Slipids_2020.ff'):
    """
    
    Parameters
    ----------
    step5_input_pdb_path : string path
        DESCRIPTION: The path to the final step of CHARMM-GUI, pdb file 
    membrane_out_path : string path
        DESCRIPTION: The path where you want to save the membrane.
        To get the membrane the following selection was used:
    membrane_select : string, gromacs selection
        DESCRIPTION: To get the membrane the following selection is used by default:
        (group System and ! group TIP3) and ! group Protein
        Therefore, this imply that your system only has:
            TIP3P, protein and the membrane, if you have ions and other molecules, you will need to change the selection properly.       
    lipid_ff_path : TYPE, optional: string path
        DESCRIPTION. The default is '/home/ale/MY_PYTHON_PACKEGES/palmiche/data/GROMACS.ff/Slipids_2020.ff'.
                Path to the lipid force field
    Returns
    -------
    None.
    """
    # Working with directories and files
    cwd = os.getcwd()
    tmpdir = tempfile.TemporaryDirectory()
    tools.cp(lipid_ff_path, tmpdir.name, r=True)
    lipid_ff_name = os.path.basename(lipid_ff_path).split('.ff')[0]
    step5_input_pdb_path = os.path.abspath(step5_input_pdb_path)
    membrane_out_path = os.path.abspath(membrane_out_path)
    
    os.chdir(tmpdir.name)

    #Create the selection
    with open("ndx.opt", "w") as opt:
        opt.write(f"\"MEMB\" {membrane_select};")
    tools.run(f"""
              echo "q" | gmx make_ndx -f {step5_input_pdb_path} -o tmp.ndx
              gmx select -s {step5_input_pdb_path} -sf ndx.opt -n tmp.ndx -on index
              """)
    
    #deleting the line _f0_t0.000 in the file
    with open("index.ndx", "r") as index:
        data = index.read()
        data = data.replace("_f0_t0.000","")
    with open("index.ndx", "w") as index:
        index.write(data)
    
    tools.run(f"""
                export GMX_MAXBACKUP=-1;         
                echo "MEMB"| gmx editconf -f {step5_input_pdb_path} -o tmpmemb.pdb -n index.ndx
                gmx pdb2gmx -f tmpmemb.pdb -o {membrane_out_path} -water none -ff {lipid_ff_name}
                """)
    os.chdir(cwd)

def protonate(mol_path, out_ext = "xyz", obabel = "openbabel.obabel"):
    """
    For this function you need openbabel
    could be installed as:
    sudo snap install openbabel
    pip install openbabel
    depending which one you have or if you installed using other method is how you call it.
    By default the call to openbable will be openbabel.obabel, but you could change this with the keyword obabel
    """
    
    if os.path.isfile(mol_path):    
        mol_name = os.path.basename(mol_path).split(".")[0]
        mol_ext = os.path.basename(mol_path).split(".")[-1]
        if mol_ext == "pdbqt":
            #No idea, but from pdbqt to xyz the H are not placed
            tmpfile = tempfile.NamedTemporaryFile(suffix='.pdb')
            tools.run(f"""
                      {obabel} {mol_path} -O {tmpfile.name} -h
                      {obabel} {tmpfile.name} -O {mol_name}.{out_ext} -h
                      """)
            tools.rm(tmpfile.name)
        else:
            tools.run(f"{obabel} {mol_path} -O {mol_name}.{out_ext} -h")
    else:
        raise FileNotFoundError(f"The file {mol_path} doesn't exist or is not accessible")

def gen_list_order(input_mol2_path, ref_mol2_path, list_mol2_out_path = "list.lis", pymol_path_cmd = "pymol", text_editor = "gedit"):
    """
    This function will open two python sessions and a text session.
    When the list had been created, then the python program will continue
    vim it is not so recommended because you will not get the help, but it will work;
    the command : "set number" will be usefull to see the number and create the list for vim.
    Parameters
    ----------
    
    input_mol2_path : string path
        DESCRIPTION. The molecule to change
    ref_mol2_path : string path
        DESCRIPTION.The reference molecule
    list_mol2_out_path : string path
        DESCRIPTION. The default is "list.lis".
        A text file with the correct oreder.
        If in the reference molecule the atoms are in the oreder A, B, C. And in the
        input the atoms are in the order B, C, A. then the list must
        have these three lines: 3, 1, 2. What means: the position
        of the number (first line, second line, and so on) represent 
        the atom order list of the reference molecule. The numbers
        on the line, represent where the reference atoms are in the 
        input structure. The file must have just the lines with the
        specification numbers.
    pymol_path_cmd : TYPE, optional
        DESCRIPTION. The default is "pymol".
    Returns
    -------
    None.
    """
    if not os.path.isfile(input_mol2_path):
        raise FileNotFoundError(f"The file {input_mol2_path} doesn't exist or is not accessible")
    if not os.path.isfile(ref_mol2_path):
        raise FileNotFoundError(f"The file {ref_mol2_path} doesn't exist or is not accessible")
    
    string = f"""
    As you see, there are two pymol sessions and a gedit sesion:
        1- The pymol session with the mol {os.path.basename(input_mol2_path)}
        is the molecule that we want to change the atom order.
        2- The pymol session with the mol {os.path.basename(ref_mol2_path)}
        is the reference molecule.
    You need to set the "label > atom identifier > ID" in pymol
    
    If in the reference molecule the atoms are in the oreder A, B, C. And in the
    input the atoms are in the order B, C, A. then the list must
    have these three lines: 3, 1, 2. What means: the position
    of the number (first line, second line, and so on) represent 
    the atom order list of the reference molecule. The numbers
    on the line, represent where the reference atoms are in the 
    input structure. The file must have just the lines with the
    specification numbers.
    
    See the following example:
    Reference molecule                      Input molecule
          5                                     9
          |                                     |
      10--1--6                              11--10--8
          |                                     |
      11--2--7                              14--5--12
          |                                     |
      12--3--8                               7--1--6
          |                                     |
      13--4--9                              4--13--3
          |                                     |
          14                                    2
    Then the list.lis will have the following lines. In parenthesis the line number:
        (1)10
        (2)5
        (3)1
        (4)13
        (5)9
        (6)8
        (7)12
        (8)6
        (9)3
        (10)11
        (11)14
        (12)7
        (13)4
        (14)2
    
    
    
    When the {text_editor} session is closed the python programme will continue, not before.
              """
    print(string)
    tools.run(f"""
              {pymol_path_cmd} {input_mol2_path}  &
              {pymol_path_cmd} {ref_mol2_path} &
              {text_editor} {list_mol2_out_path}
              """)
    
    

def order_mol2(input_mol2_path, ref_mol2_path, list_mol2_path, lname_mol2 = 'UNL', obabel = "openbabel.obabel"):

    """
    
    Parameters
    ----------
    input_mol2_path : string, path
        The molecule that you want to change the atom list
        order in any openbabel supported format. Is important that 
        it has the same number of atoms that the reference. Check
        for the H atoms, if the implicit H are missing, run first
        obabel mol.abc -O mol.xyz -h
        and compare the structures in Pymol in order to 
        create the list file.
        
    ref_mol2_path : string, path
        The molecule that you will use as reference for the 
        atom list order in mol2 format
    
    list_mol2_path : string, path
        A text file with the correct oreder. If in the ref
        structure the atoms are in the oreder A, B, C. And in the
        input the atoms are in the order B, C, A. then the list must
        have these three lines: 3, 1, 2. What means: the position
        of the number (first line, second line, and so on) represent 
        the atom order list of the reference structure. The numbers
        on the line, represents where the reference atom is in the 
        input structure. The file must have just the lines with the
        specification   numbers.
    
    lname_mol2 : string, optional
       The three letter identifier for the ligand. Default is UNL.
    obabel : string, optional
       The name for execute openbabel in case that the extension of input_mol2_path is different from .xyz . Default is openbabel.obabel.
    Returns
    -------
    None.
    """
    file_input_mol2 = os.path.basename(input_mol2_path)
    out_path = os.path.dirname(input_mol2_path)
    
    #print(os.path.exists(input_mol2_path), input_mol2_path)
    if file_input_mol2.split('.')[-1] == 'xyz':
        to_work = input_mol2_path
    else:
        xyztmp = tempfile.NamedTemporaryFile(suffix=".xyz")
        tools.run(f'{obabel} {input_mol2_path} -O {xyztmp.name}')
        to_work = xyztmp.name
    
    order = [int(i) - 1 for i in open(list_mol2_path, 'r').readlines()] # -1 is because the python array start at 0
    
    with open(to_work, 'r') as f:
        lines_input = f.readlines()
    Natoms = int(lines_input[0])
    coords = lines_input[2:2+Natoms]

    if len(coords) != len(order): exit(f"the len of coords in {input_mol2_path} doesn't match with the number of line in {ref_mol2_path}")
    
    update_coords = []
    for i in order:
        update_coords.append(coords[i])
    update_coords = [c.split()[1:] for c in update_coords]
    
    with open(ref_mol2_path, 'r') as f:
        lines_ref = f.readlines()
    name_final_out = file_input_mol2.split('.')[0]+".mol2"
    
    #The file is written to the directory where the input is
    with open(os.path.join(out_path, name_final_out), 'wt') as final_out:
        i = 0
        while i < len(lines_ref):
            final_out.write(lines_ref[i])
            if lines_ref[i].startswith('@<TRIPOS>MOLECULE'):
                final_out.write(f'    {lname_mol2}\n')
                i += 1
            if lines_ref[i].startswith('@<TRIPOS>ATOM'):
                for j in range(len(update_coords)):
                    
                    atom_split = lines_ref[i+j+1].split()
                    to_print = "\t".join(atom_split[:2]) + "\t" + "\t".join(update_coords[j]) + '\t' + "\t".join(atom_split[-4:-2]) + "\t" + lname_mol2 + "\t" + atom_split[-1] + '\n'   
                    
                    final_out.write(to_print)
                i = i + Natoms
            i += 1  

def run_tleap(mol2_path, frcmod_path):
    name = os.path.basename(mol2_path).split(".")[0]    
    tmpin = tempfile.NamedTemporaryFile(suffix='.in')

    with open(tmpin.name, 'w') as t:
        t.write(f"""
                  loadAmberParams {frcmod_path}
                  A = loadMol2 {mol2_path}
                  saveAmberParm A {name}.prmtop {name}.crd
                  quit
                  """)
        
    tools.run(f"""
               tleap -s -f {tmpin.name} > {name}_tleap.out
               acpype -p {name}.prmtop -x {name}.crd
               """)
    
    [tools.run(f"obabel {gro} -O {gro.split('.')[0]}.pdb") for gro in glob("*_GMX.gro")]



def toptleap2itp(topology, out_prefix = "UNL"):
    sections = tools.get_top_sections(topology)
    
    #Getting the important sections, in a dict is much more easy, but we
    #have the problem of the duplication of dihedrals section
    for section in sections:
        if section[0] == "atoms":
            atom_section = section
        elif section[0] == "atomtypes":
            atomtypes_section = section
        else:
            pass
    
    #I need to run again the loop because I am not able to pop or delete on 
    #the fly, because the for loop will raise error, I am not able to sort
    #Because I dont know the implications of it.
    itp_sections = ""
    for section in sections:
        if section[0] not in ["defaults" , "atomtypes", "system", "molecules"]:
            itp_sections += f"[ {section[0]} ]\n{''.join(section[1])}"
        

    #Writing the atomtypes itp
    with open("atomtypes.itp", "w") as at:
        at.write(f"[ {atomtypes_section[0]} ]\n{''.join(atomtypes_section[1])}")
        
    #Writing the general itp
    #Even when the five ligand on different chains will have the same
    #contents on this file. You need to make a distinction because the
    #moleculetype definition inside of the topology file that references
    #the name of the molecule. And this is needed for the pulling simulations
    
    include_posre = f"""
; Include Position restraint file for heavy atoms
#ifdef POSRES
#include "{out_prefix}_posre.itp"
#endif
"""
    with open(f"{out_prefix}.itp", "w") as itp:
        itp.write(f"{itp_sections}{include_posre}")
    
    #creating the posre, first we will look for the heavy atoms
    heavy_atoms = []
    for atom in atom_section[1]:
        if not atom.startswith(";"):#This is necessary because the commented lines in gromacs topology files
            split = atom.strip().split()#strip for the case of blank lines
            if split:#The list has something, this is all for the case of blank lines
                if split[4][0].lower() != "h": #This is the column that has the atom name, I select the first char because could be H0, HH13, etc, and all of them are Hydrogens atoms, hoping that all of them start with H
                    heavy_atoms.append(split[0])#Here is the column that has the atom index
            
    #Writing the posre file
    with open(f"{out_prefix}_posre.itp", "w") as posre:
        posre.write("[ position_restraints ]\n")
        posre.write(f"{';i':<5}{'funct':^5}{'fcx':^20}{'fcy':^20}{'fcz':^20}\n")
        for atom_index in heavy_atoms:
            posre.write(f"{atom_index:^5}{1:^5}{'POSRES_LIG':^20}{'POSRES_LIG':^20}{'POSRES_LIG':^20}\n")
                
        
def topoltop(topol, numb_lipids, out_name = None, backup = True):
    """
    !!!!!This function is NOT general at all, it is just created here for stetic
    and is strong depended of what I needed to do in my project, so, use it
    with precaution
    Parameters
    ----------
    topol : TYPE, path, string
        DESCRIPTION. a path or name to a GROMACS topology file
    out_name : TYPE, string
        DESCRIPTION. The name of the output, if None, teh topol will be used
        It is just the name, not the path, the path will be the path were the topol file is
    Returns
    -------
    None.
    the modified topology
    """               
    ff_include_old = """
; Include forcefield parameters
#include "./amber99sb-star-ildn.ff/forcefield.itp"
"""
    ff_include_new = """
; Include forcefield parameters
#include "./amber99sb-star-ildn.ff/forcefield.itp"
; Include forcefield parameters		
#include "./Slipids_2020.ff/forcefield.itp"
; Include GAFF2 atom definition
#include "./GAFF2.ff/atomtypes.itp"
"""

    rest_str_old = """
; Include water topology
#include "./amber99sb-star-ildn.ff/tip3p.itp"
#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif
; Include topology for ions
#include "./amber99sb-star-ildn.ff/ions.itp"
[ system ]
; Name
Protein
[ molecules ]
; Compound        #mols
Protein_chain_A     1
Protein_chain_B     1
Protein_chain_C     1
Protein_chain_D     1
Protein_chain_E     1"""

    rest_str_new = f"""
; Include POPC chain topology
#include "toppar/POPC.itp"
; position restraints for POPC lipid
#ifdef POSRES
[ position_restraints ]
   20     1     0.0             0.0            POSRES_FC_LIPID
#endif
#ifdef DIHRES
[ dihedral_restraints ]
   25    36    28    30     1   -120.0      2.5       DIHRES_FC
   60    63    65    67     1      0.0      0.0       DIHRES_FC
#endif	
; Include ligand topologies
#include "toppar/LIA.itp"
#include "toppar/LIB.itp"
#include "toppar/LIC.itp"
#include "toppar/LID.itp"
#include "toppar/LIE.itp"
; Include water topology
#include "./amber99sb-star-ildn.ff/tip3p.itp"
#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif
; Include topology for ions
#include "./amber99sb-star-ildn.ff/ions.itp"
[ system ]
; Name
Protein
[ molecules ]
; Compound        #mols
Protein_chain_A     1
Protein_chain_B     1
Protein_chain_C     1
Protein_chain_D     1
Protein_chain_E     1
POPC		{numb_lipids}
LIA		1
LIB		1
LIC		1
LID		1
LIE		1"""
   
    with open(topol, "r") as t:
        topol_data = t.read()
    topol_data = topol_data.replace(ff_include_old, ff_include_new)
    topol_data = topol_data.replace(rest_str_old, rest_str_new)
    topol_data = topol_data.replace("topol_Protein_chain_", "toppar/topol_Protein_chain_")
       
    if out_name:
        with open(os.path.join(os.path.dirname(topol), out_name), "w") as out:
            out.write(topol_data)
    else:
        if backup: tools.backoff(topol)
        with open(topol, "w") as out:
            out.write(topol_data)

def posre_BB_SC(pdb_path):
    tools.run(f"""
    subs='1e+08'
    
    echo "Backbone" | gmx genrestr -f {pdb_path} -o backbone_posre -fc $subs $subs $subs
    sed "s/$subs/POSRES_FC_BB/g" backbone_posre.itp > tmp ; mv tmp backbone_posre.itp
    
    echo "SideChain-H" | gmx genrestr -f {pdb_path} -o sidechain_posre -fc $subs $subs $subs
    sed "s/$subs/POSRES_FC_SC/g" sidechain_posre.itp > tmp ; mv tmp sidechain_posre.itp
              """)
    
def include_posre_itp(itp_in, *args, itp_out = None, backup =True):
    """
    
    Parameters
    ----------
    itp_in : string path
        DESCRIPTION: .itp file form GROMACS
    *args : string
        DESCRIPTION: You must give at least one,
        These are the file that you will included with the position restrain.
        See the example below:
           
            ; Include Position restraint file for protein backbone
            #ifdef POSRES
            #include "{arg}"
            #endif
            
    itp_out : TYPE, optional, string path
        DESCRIPTION. The default is None. The path to the output file, if None,
        the original file will be updated and the former backed it up
    Raises
    ------
    KeyError
        DESCRIPTION. In case that you don't give extra arguments
    Returns
    -------
    None.
    """
    if not len(args):
        raise KeyError("You need to provied at least on argument to fill in the itp")
    
    include_str = ""
    for arg in args:
        include_str += f"""
; Include Position restraint file for protein
#ifdef POSRES
#include "{arg}"
#endif
        """
     
    with open(itp_in, "r") as f:
        itp_lines = f.readlines()
    for (i,line) in enumerate(reversed(itp_lines)):
        if "; Include Position restraint file" in line:
            break
    #Here in this line are bassically eliminating the position restraing definition
    #And replaced by a new one
    itp_lines = "".join(itp_lines[:len(itp_lines) - 1 - i]) + include_str
    
    if itp_out:
        with open(itp_out, "w") as f:
            f.write(itp_lines)
    else:
        if backup: tools.backoff(itp_in)
        with open(itp_in, "w") as f:
            f.write(itp_lines)       
def gmxsolvate(abs_path, input_file, topology, editconf_box, editconf_angles = (90,90,60), editconf_bt = "tric", solvate_cs = "spc216", out_file = "solvated.pdb"):
    """
    
    Parameters
    ----------
    abs_path : TYPE, string path
        DESCRIPTION: Here we must deffine the path were all the necessary files
        are together, topologies, forcefields, configurations files, etc
    input_file : TYPE, string
        DESCRIPTION. Name of the configuration input file
    topology : TYPE, string
        DESCRIPTION. GROMACS topology file
    editconf_box : TYPE tuple
        DESCRIPTION. HEre we deffine the components of the vector on editconf (a,b,c)
    editconf_angles : TYPE, optional, tuple
        DESCRIPTION. The default is (90,90,60). This default value is for
        hexagonal boxes. idela for membrane simulations
    editconf_bt : TYPE, optional, string
        DESCRIPTION. The default is "tric". Type of GROMACS box
    solvate_cs : TYPE, optional, string
        DESCRIPTION. The default is "spc216". Configuration for the solvent
    out_file : TYPE, optional, string
        DESCRIPTION. The default is "solvated.pdb". Name of the solvated box
    Returns
    -------
    None.
    """
    
    initial_cwd = os.getcwd()
    
    os.chdir(abs_path)
    tools.run(f"""
                export GMX_MAXBACKUP=-1
                gmx editconf -f {input_file} -o {out_file} -bt {editconf_bt} -box {' '.join([str(i) for i in editconf_box])} -angles {' '.join([str(i) for i in editconf_angles])}
                gmx solvate -cp {out_file} -p {topology} -cs {solvate_cs} -o {out_file}
              """)
             
    os.chdir(initial_cwd)

def gmxions(abs_path, input_file, topology, out_file = "ions.pdb", genion_pname = "K", genion_nname = "CL", genion_rmin = 1.0):
    """
    
    Parameters
    ----------
    abs_path : TYPE
        DESCRIPTION.
    input_file : TYPE, string 
        DESCRIPTION. the name of the input file inside the abs_path
    topology : TYPE, string
        DESCRIPTION. the name of the topology file inside the abs_path
    out_file : TYPE, optional, string
        DESCRIPTION. The default is "ions.pdb". The name of the output file
        inside the abs_path. Path will not handlead.
        
        For the rest of the parameters see the fromacs documentation of genion
    genion_pname : TYPE, optional
        DESCRIPTION. The default is "K".
    genion_nname : TYPE, optional
        DESCRIPTION. The default is "CL".
    genion_rmin : TYPE, optional
        DESCRIPTION. The default is 1.0.
    Returns
    -------
    None.
    """
    initial_cwd = os.getcwd()
    tmptpr = tempfile.NamedTemporaryFile(suffix = '.tpr')
    tmpmdp = tempfile.NamedTemporaryFile(suffix = '.mdp')
    os.chdir(abs_path)
    tools.run(f"""
                export GMX_MAXBACKUP=-1
                gmx grompp -f {tmpmdp.name} -c {input_file} -p {topology} -o {tmptpr.name}
                echo "SOL" | gmx genion -s {tmptpr.name} -p {topology} -o {out_file} -neutral -pname {genion_pname} -nname {genion_nname} -rmin {genion_rmin}
              """)
    # Cleanning
    tools.rm("mdout.mdp")
    os.chdir(initial_cwd)

def finalndx(abs_path, input_file, ndxout = "index.ndx", chainsID = ["A", "B", "C", "D", "E"]):
    
    initial_cwd = os.getcwd()
    os.chdir(abs_path)
    
    tmpopt = tempfile.NamedTemporaryFile(suffix='.opt')
    tmpndx = tempfile.NamedTemporaryFile(suffix='.ndx')
    #Nice use of gmx select, see the use of the parenthesis
    with open(tmpopt.name, "w") as opt:
        opt.write(f"""
                  "MEMB" ((group System and ! group Water_and_ions) and ! group Protein) and ! ({"resname LI"+" or resname LI".join(chainsID)});
                  "SOLU" group Protein or {"resname LI"+" or resname LI".join(chainsID)};
                  "SOLV" group Water_and_ions;
                  """)
    tools.run(f"""
                export GMX_MAXBACKUP=-1
                echo "q" | gmx make_ndx -f {input_file} -o {tmpndx.name}
                gmx select -s {input_file} -sf {tmpopt.name} -n {tmpndx.name} -on {ndxout}
                """)
    
    #deleting the line _f0_t0.000 in the file
    with open(ndxout, "r") as index:
        data = index.read()
        data = data.replace("_f0_t0.000","")
    with open(ndxout, "w") as index:
        index.write(data)
    
    
    os.chdir(initial_cwd)

def GROMACS_SMAUG_JOB(job_name, conf_file_name = "ASSEMBLE.pdb", 
                               partition="deflt", 
                               output="myjob.out", error="myjob.err",
                               cpus_per_task=12, gpus = 0, time="2-00:00",
                               nodes=1, nice=0, version = "2021.1", exclude=None):
    template = f"""#!/bin/bash
#SBATCH --partition {partition}
#SBATCH --output {output}
#SBATCH --error {error}
#SBATCH --cpus-per-task={cpus_per_task}
#SBATCH --time {time}
#SBATCH --job-name={job_name}
#SBATCH --nodes={nodes}
#SBATCH --nice={nice}
#SBATCH --gpus={gpus}
    """
    if exclude:
        template += f"#SBATCH --exclude={exclude}\n"
    template += f"""
\n\n
# This block is echoing some SLURM variables
echo "Job execution start: $(date)"
echo "JobID = $SLURM_JOBID"
echo "Host = $SLURM_JOB_NODELIST"
echo "Jobname = $SLURM_JOB_NAME"
echo "Subcwd = $SLURM_SUBMIT_DIR"
echo "SLURM_TASKS_PER_NODE = $SLURM_TASKS_PER_NODE"
echo "SLURM_CPUS_PER_TASK = $SLURM_CPUS_PER_TASK"
echo "SLURM_CPUS_ON_NODE = $SLURM_CPUS_ON_NODE"
cd $(pwd)
source /data/shared/opt/gromacs/{version}/bin/GMXRC
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:
# Minimization
gmx grompp -f step6.0_minimization.mdp -o step6.0_minimization.tpr -c {conf_file_name} -r {conf_file_name} -p topol.top -n index.ndx
gmx mdrun -v -deffnm step6.0_minimization
# Equilibration
cnt=1
cntmax=6
for ((i=${{cnt}}; i<${{cntmax}}+1; i++)); do
	if [ $i == "1" ]; then
		gmx grompp -f step6.${{i}}_equilibration.mdp -o step6.${{i}}_equilibration.tpr -c step6.$[${{i}}-1]_minimization.gro -r {conf_file_name} -p topol.top -n index.ndx
		gmx mdrun -v -deffnm step6.${{i}}_equilibration -nt 12
	else
		gmx grompp -f step6.${{i}}_equilibration.mdp -o step6.${{i}}_equilibration.tpr -c step6.$[${{i}}-1]_equilibration.gro -r {conf_file_name} -p topol.top -n index.ndx
		gmx mdrun -v -deffnm step6.${{i}}_equilibration -nt 12
	fi
done 
    """
    return(template)

if __name__ == "__main__":
    pass

            
