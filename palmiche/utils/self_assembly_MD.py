#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from palmiche.utils import tools, jobsh
from palmiche.utils.tools import PathLike
from palmiche.utils import pdb, top
from palmiche.utils.pKa2pdb2gmx import pka2gmx
from glob import glob
import tempfile
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
import openbabel
from toff import Parameterize


def get_path_dict(ligands_path:PathLike = "./ligands", receptors_path:PathLike = "./receptors", dockings_path:PathLike = "./run_vina_run"):

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
    # for directory in tools.list_if_dir(ligands_path):
    for file in tools.list_if_file(ligands_path):
        file_name, _ = os.path.splitext(file)
        if os.path.isdir(os.path.join(ligands_path, file_name)):
            ligand_dict[file_name] = {"frcmod": os.path.join(ligands_path, file_name, "mol.frcmod"),
                                    "mol2": os.path.join(ligands_path, file_name, "mol.mol2"),
                                    "pdb": os.path.join(ligands_path, file_name, "mol.pdb"),
                                    "list": os.path.join(ligands_path, file_name, "list.lis")}
        else:
            ligand_dict[file_name] = None

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

def fix_pdb(step5_input_pdb:PathLike, out:PathLike = "PROTEIN.pdb"):
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

    cwd = os.getcwd()
    tmp_dir = tempfile.TemporaryDirectory(prefix='.fix_pdb_', dir = cwd)
    os.chdir(tmp_dir.name)

    _, file_extension = os.path.splitext(step5_input_pdb)
    if file_extension == "gro":
        tools.run(f"export GMX_MAXBACKUP=-1; echo gmx editconf -f {step5_input_pdb} -o tmp1.pdb")

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
        tools.cp(step5_input_pdb, 'tmp1.pdb')

    #Keeping only the protein
    tools.run(f"""
            export GMX_MAXBACKUP=-1
            echo "q" | gmx make_ndx -f tmp1.pdb -o tmp.ndx
            echo "Protein" | gmx editconf -f tmp1.pdb -o tmp2.pdb -n tmp.ndx
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


    PDB_OBJECT = pdb.PDB('tmp2.pdb')
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
    PDB_OBJECT.write('tmp2.pdb', backup = False)
    #Get correct protonation, here I use amber99sb-ildn because it is by default in GROMACS

    #Probably a better idea is use CHARMM force fields in irder to dont worry about the atom names, but this use (amber) is consistent with my project
    pka2gmx('tmp2.pdb',
            output = 'tmp2.pdb',
            protein_forcefield = "amber99sb-ildn",
            water_forcefield = "tip3p",
            ph = 7.0,
            pKa = "propka3",
            pKa_file_out = False,
            itp_top_out = False,
            num_file = False)
    PDB_OBJECT = pdb.PDB('tmp2.pdb') #Need to load the new pdb

    PDB_OBJECT.write(out, backup = False)
    os.chdir(cwd)
    tmp_dir.cleanup()
    return PDB_OBJECT

def get_cryst1(step5_input_gro_path):
    with tempfile.NamedTemporaryFile(suffix='.pdb') as tmp_pdb:
        tools.run(f"gmx editconf -f {step5_input_gro_path} -o {tmp_pdb.name}")
        TMP_PDB = pdb.PDB(tmp_pdb.name)
        #I need to dived by 10 because the pdb is in Angstrom
        box_vector = (TMP_PDB.cryst1.a/10, TMP_PDB.cryst1.b/10, TMP_PDB.cryst1.c/10)
        box_angles = (TMP_PDB.cryst1.alpha, TMP_PDB.cryst1.beta, TMP_PDB.cryst1.gamma)
    return box_vector, box_angles


def gmxmemb(step5_input_pdb_path,
            membrane_out_path,
            membrane_select = "(group System and ! group TIP3) and ! group Protein",
            lipid_ff_path = 'Slipids_2020.ff'):
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
        DESCRIPTION. The default is 'Slipids_2020.ff'.
                Path to the lipid force field
    Returns
    -------
    None.

    """
    # Working with directories and files
    cwd = os.getcwd()
    tmpdir = tempfile.TemporaryDirectory(prefix='.gmxmemb_', dir=cwd)
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
    tmpdir.cleanup()



def run_tleap(mol2_path, frcmod_path):
    file_name, _ = os.path.splitext(mol2_path)
    with tempfile.NamedTemporaryFile(suffix='.in') as tmpin:
        with open(tmpin.name, 'w') as t:
            string = f"LoadAmberParams {frcmod_path}\n"\
                    f"A = loadMol2 {mol2_path}\n"\
                    f"saveAmberParm A {file_name}.prmtop {file_name}.crd\n"\
                    "quit"
            t.write(string)

        tools.run(f"""
                tleap -s -f {tmpin.name} > {file_name}_tleap.out
                acpype -p {file_name}.prmtop -x {file_name}.crd
                """)


# TODO A much general function here. the problem is in how to read the include statements
def topoltop(topol, numb_lipids, out_name = None, ligand_ff_code = 'GAFF2', backup = True):
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
    ff_include_new = f"""
; Include forcefield parameters
#include "./amber99sb-star-ildn.ff/forcefield.itp"

; Include forcefield parameters
#include "./Slipids_2020.ff/forcefield.itp"

; Include {ligand_ff_code} atom definition
#include "./{ligand_ff_code}.ff/atomtypes.itp"
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
    tmptpr.close()
    tmpmdp.close()
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

    tmpopt.close()
    tmpndx.close()
    os.chdir(initial_cwd)

def GROMACS_JOB(output = 'job.sh', conf='ASSEMBLY.pdb', hostname = 'smaug', GROMACS_version = '2021.5'):
    GROMACS_section = f"""
# Minimization

gmx grompp -f step6.0_minimization.mdp -o step6.0_minimization.tpr -c {conf} -r {conf} -p topol.top -n index.ndx
gmx mdrun -v -deffnm step6.0_minimization

# Equilibration
cnt=1
cntmax=6

for ((i=${{cnt}}; i<${{cntmax}}+1; i++)); do
	if [ $i == "1" ]; then
		gmx grompp -f step6.${{i}}_equilibration.mdp -o step6.${{i}}_equilibration.tpr -c step6.$[${{i}}-1]_minimization.gro -r {conf} -p topol.top -n index.ndx
		gmx mdrun -v -deffnm step6.${{i}}_equilibration -nt 12
	else
		gmx grompp -f step6.${{i}}_equilibration.mdp -o step6.${{i}}_equilibration.tpr -c step6.$[${{i}}-1]_equilibration.gro -r {conf} -p topol.top -n index.ndx
		gmx mdrun -v -deffnm step6.${{i}}_equilibration -nt 12
	fi
done
    """
    jobsh_obj = jobsh.JOB(
    sbatch_keywords={'job-name':'assamble'},
    GROMACS_version=GROMACS_version,
    build_GROMACS_section=GROMACS_section,
    hostname=hostname,
    )
    jobsh_obj.write(output)

def MolFromPdbqtFile(pdbqt:PathLike) -> Chem.rdchem.Mol:
    """Return a protonanted RDKit molecule.

    Parameters
    ----------
    pdbqt : PathLike
        The path to the pdbqt file.

    Returns
    -------
    Chem.rdchem.Mol
        The molecule.
    """
    pdbqt_mol = PDBQTMolecule.from_file(pdbqt, skip_typing=True, poses_to_read=1)
    mol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0]
    if not mol:
        pdb_tmp = tempfile.NamedTemporaryFile(suffix='.pdb')

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("pdbqt", "pdb")
        mol_obabel = openbabel.OBMol()
        obConversion.ReadFile(mol_obabel, pdbqt)
        obConversion.WriteFile(mol_obabel, pdb_tmp.name)

        # I am not sure why not ony the following line
        # For some reason Parameterize fails, I think that may be some initialization of Molecule
        # TODO check here
        # It generates in TOFF
        # Exception: Topology and System have different numbers of atoms (10 vs. 24)
        # It is something releated with the pdb reader
        # mol = Chem.MolFromPDBFile(pdb_tmp.name)
        mol0 = Chem.MolFromPDBFile(pdb_tmp.name)
        mol0 = Chem.AddHs(mol0, addCoords = True)
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol0))
        mol = Chem.AddHs(mol, addCoords = True)
        tools.replace_conformer(mol, mol0)

        pdb_tmp.close()
    mol = Chem.AddHs(mol, addCoords = True)
    Chem.AssignStereochemistryFrom3D(mol)
    # I am not sure Why I need this, but the molecule readed from the dp
    # to_return = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    return mol

def make_ligand_topology(
    receptor_code:str,
    ligand_code:str,
    docking_dict:dict,
    ligand_dict:dict,
    ligand_ff_code:str = 'GAFF2',
    openff_code:str = 'openff_unconstrained-2.0.0.offxml',
    out_dir:PathLike = '.'):

    #First, create the directory
    cwd = os.getcwd()
    tmp_dir = tempfile.TemporaryDirectory(prefix=f'.make_lig_top_{receptor_code}_{ligand_code}_{ligand_ff_code}', dir=cwd)
    os.chdir(tmp_dir.name)

    ligand_dir = os.path.join(out_dir, receptor_code, ligand_code)
    tools.makedirs(ligand_dir)

    toppar_dir = os.path.join(ligand_dir, "toppar")
    tools.makedirs(os.path.join(toppar_dir))

    ligand_ff_dir = (os.path.join(ligand_dir, f"{ligand_ff_code}.ff"))
    tools.makedirs(ligand_ff_dir)

    #Here we are selected the BE pose, because the results were really really god and the orientation was the same. A perfect result
    #With sorted we ensure the order of the chains
    tools.touch(os.path.join(ligand_dir, "ASSEMBLY.pdb"))#Create first an empty file
    ASSEMBLE_PDB = pdb.PDB(os.path.join(ligand_dir, "ASSEMBLY.pdb"))#Empty pdb object, take the former empy file
    for i, chain in enumerate(sorted(docking_dict[receptor_code][ligand_code])):
        #Here we are taking the Best Energy conformation
        name = os.path.basename(docking_dict[receptor_code][ligand_code][chain]["BE"]).split(".")[0]
        # Get the docking conformation
        docking_mol = MolFromPdbqtFile(docking_dict[receptor_code][ligand_code][chain]["BE"])
        if ligand_ff_code == 'GAFF2':
            # Add Docking conformation and change name to mol2
            mol2 = tools.Mol2(ligand_dict[ligand_code]["mol2"])
            mol2.AddConformer(docking_mol)
            mol2.ChangeName(f"LI{chain}")
            mol2.write(f"{name}.mol2")
            # Create the topology
            run_tleap(f"{name}.mol2",ligand_dict[ligand_code]["frcmod"])

            # Convert gro to pdb
            gro_file = f"LI{chain}.amb2gmx/LI{chain}_GMX.gro"
            obConversion = openbabel.OBConversion()
            obConversion.SetInAndOutFormats("gro", "pdb")
            mol_obabel = openbabel.OBMol()
            obConversion.ReadFile(mol_obabel, gro_file)
            obConversion.WriteFile(mol_obabel,f"LI{chain}.pdb")

            # Get topology
            tools.mv(os.path.join(f"LI{chain}.amb2gmx", f"LI{chain}_GMX.top"),f"./LI{chain}.top")
            # Here I have to add the condition when OpenFF is used
        elif ligand_ff_code == 'OpenFF':
            paramaterizer = Parameterize(
                force_field_code = openff_code,
                ext_types = ['top', 'pdb'],
                hmr_factor = None,
                overwrite = True,
                safe_naming_prefix = 'z',
                out_dir = '.',
            )
            paramaterizer(input_mol=docking_mol, mol_resi_name=f"LI{chain}")
            # print(os.listdir('.'))
            # raise RuntimeError


        #Here we are concatenating all the atom objects of each ligand in one

        tmp_LIGAND = pdb.PDB(f"LI{chain}.pdb")
        #This is nedded because if not all the ligands have the same chain "A"
        for atom in tmp_LIGAND.atoms:
            atom.chainID = chain

        ASSEMBLE_PDB.atoms += tmp_LIGAND.atoms

        # Convert the topology to itp
        ligand_topology = top.TOP(f"LI{chain}.top")
        ligand_topology.to_itp()
        ligand_topology.make_posres(posre_name = f"LI{chain}", out_dir=toppar_dir)

        ligand_topology.write(os.path.join(toppar_dir, f"LI{chain}.itp"))

        # Get the atomtypes section only one time
        if i == 0:
            ligand_topology._read()
            for section in list(ligand_topology.get_section_names()):
                if section != 'atomtypes':
                    try:
                        ligand_topology.del_section(section)
                    except KeyError:
                        # this esception reise in case of repeated sections
                        pass
            ligand_topology.write(os.path.join(ligand_ff_dir,"atomtypes.itp"))

        # [tools.mv(item, toppar_dir) for item in [f"LI{chain}.itp", f"LI{chain}_posre.itp"]]

    os.chdir(cwd)
    tmp_dir.cleanup()
    return ASSEMBLE_PDB
    #After we create all the ligands per chain we write the result in the
    #corresponding directory but also with the atoms of protein, membrane + ligands
    #Becasue we dont especified, it will use the directory used for its construction


def assembly(
    ligands_path_root:PathLike = './ligand',
    receptors_path_root:PathLike = './receptor',
    dockings_path_root:PathLike = './run_vina_run',
    receptor_codes:list[str] = None,
    ligand_codes:list[str] = None,
    ligand_ff_code = 'GAFF2', # OpenFF
    openff_code:str = 'openff_unconstrained-2.0.0.offxml',
    out_dir = 'Assembly',
    hostname = 'smaug',
    GROMACS_version = '2021.5',
    increase_box_on_z = None,

    ):

    # Create a temporal working directory
    tmp_wd = tempfile.TemporaryDirectory(prefix='.selfassembly', dir='.')
    # Convert to absolute path
    ligands_path_root = os.path.abspath(ligands_path_root)
    receptors_path_root = os.path.abspath(receptors_path_root)
    dockings_path_root = os.path.abspath(dockings_path_root)
    out_dir = os.path.abspath(out_dir)
    tools.makedirs(out_dir)
    cwd = os.getcwd()
    os.chdir(tmp_wd.name)

    ff_code_translate = {
        'gaff2': 'GAFF2',
        'openff': 'OpenFF',
    }
    ligand_ff_code = ff_code_translate[ligand_ff_code.lower()]
    if ligand_ff_code not in ff_code_translate.values():
        raise RuntimeError(f"Non valid ligand_ff_code. Choose gaff2 or openff (case insensitive).")
    #Deffinitions of the force fields
    water_ff = "tip3p"
    propka = 'propka3'


    receptor_dict, ligand_dict, docking_dict = get_path_dict(
        ligands_path = ligands_path_root,
        receptors_path = receptors_path_root,
        dockings_path = dockings_path_root
    )
    if receptor_codes:
        for receptor in list(receptor_dict):
            if receptor not in receptor_codes:
                del receptor_dict[receptor]
                if receptor in docking_dict:
                    del docking_dict[receptor]
    if ligand_codes:
        for receptor in list(docking_dict):
            for ligand in list(docking_dict[receptor]):
                if ligand not in ligand_codes:
                    del docking_dict[receptor][ligand]
                    del ligand_dict[ligand]

    for receptor in receptor_dict:

        #Adding new keys to receptor_dict
        dirname = os.path.dirname(receptor_dict[receptor]['protein'])
        receptor_dict[receptor]['update'] = os.path.join(dirname,f'{receptor}_update.pdb')
        receptor_dict[receptor]['topol'] = os.path.join(dirname, 'topol.top')
        receptor_dict[receptor]['pKa'] = os.path.join(dirname, f'{receptor}.pka')
        receptor_dict[receptor]['toppar'] = os.path.join(dirname, 'toppar/')

        #!!Get the vector information in the step5_input.gro of CHARMM-GUI, this info is only in the gro file
        vector, angles = get_cryst1(receptor_dict[receptor]['step5_input_gro'])
        if increase_box_on_z:
            vector = list(vector)
            vector[2] += increase_box_on_z
            
        receptor_dict[receptor]['cryst1'] = {'vector':vector, 'angles':angles}

        #If the membrane is not already created (looking the path storged in
        #protein["membrane"]), processed with pdb2gmx
        if not os.path.isfile(receptor_dict[receptor]['membrane']):
            gmxmemb(receptor_dict[receptor]['step5_input_pdb'], receptor_dict[receptor]['membrane'], lipid_ff_path = lipid_ff_path)


        #This alwayas should be ignored if you came from the docking steps
        if not os.path.isfile(receptor_dict[receptor]['protein']):
            raise FileNotFoundError(f"The file {receptor_dict[receptor]['protein']} doesn't exist, you should construct the protein aligned to the docking results. Use vina_auto first.")

        #Storage the PROT_PDB and chains
        PROT_PDB = pdb.PDB(receptor_dict[receptor]['protein'])

        #Update the protein with the corresponded force field
        #Probably I also need to consider use the final force field on that step to
        #avoid this second evaluation of pka2gmx
        #Only if the protein_update.pdb doesn't exist
        if not os.path.isfile(receptor_dict[receptor]["update"]) or not os.path.isfile(receptor_dict[receptor]["topol"]) or not os.path.isdir(receptor_dict[receptor]["toppar"]):
            tools.get_palmiche_data(file="GROMACS.ff/amber99sb-star-ildn.ff.tar.gz", out_dir='.')
            tools.get_palmiche_data(file="GROMACS.ff/Slipids_2020.ff.tar.gz", out_dir='.')
            tools.get_palmiche_data(file="MDP/CHARMM-GUI", out_dir='.')
            pka2gmx(receptor_dict[receptor]["protein"],
                    output = f"{receptor}_update.pdb",
                    protein_forcefield = 'amber99sb-star-ildn.ff',
                    water_forcefield = water_ff,
                    ph = 7.0,
                    pKa = propka)
            tools.mv(f"{receptor}.pka", receptor_dict[receptor]["pKa"])
            tools.mv(f"{receptor}_update.pdb", receptor_dict[receptor]["update"])



            #!!!Modifify the topo.top. The fucntion topoltop IS NOT general and is
            #strong depended of the system, force, field, etc.... Take a look of the
            #function and WORKX WITH PRECAUTION!.
            #make a save of the roiginal topol.top

            tools.cp("topol.top", os.path.join(dirname, f"original_topol_{receptor}_update.top"))
            numb_lipids = len(pdb.PDB(receptor_dict[receptor]["membrane"]).get_residues())
            topoltop("topol.top", numb_lipids, ligand_ff_code = ligand_ff_code, out_name = None)
            tools.mv("topol.top", receptor_dict[receptor]["topol"])

            #Create the side chain and back bone posre itp files
            for chain in PROT_PDB.get_chains():
                tmp_pdb = 'tmp_posre.pdb'
                tools.touch(tmp_pdb)
                TMP_PDB = pdb.PDB(tmp_pdb)#Empty pdb object, take the former empty file
                TMP_PDB.atoms = [atom for atom in PROT_PDB.atoms if atom.chainID == chain]
                TMP_PDB.write(backup=False)
                posre_BB_SC(tmp_pdb)
                os.rename("backbone_posre.itp", f"backbone_posre_{chain}.itp")
                os.rename("sidechain_posre.itp", f"sidechain_posre_{chain}.itp")
                #Include in all the itp the corresponding posre files
                include_posre_itp(f"topol_Protein_chain_{chain}.itp", f"backbone_posre_{chain}.itp", f"sidechain_posre_{chain}.itp")
            #Storoge and create a path to the itp files
            tools.makedirs(receptor_dict[receptor]["toppar"])
            [tools.mv(item, receptor_dict[receptor]["toppar"]) for item in glob(f"topol_Protein_chain_{PROT_PDB.get_chains()}.itp")+glob(f"backbone_posre_{PROT_PDB.get_chains()}.itp")+glob(f"sidechain_posre_{PROT_PDB.get_chains()}.itp")]

        # Change the code of the of the force field in the topology
        if ligand_ff_code == 'OpenFF':
            possible_old_ff_code = 'GAFF2'
        else:
            possible_old_ff_code = 'OpenFF'
        with open(receptor_dict[receptor]["topol"], 'r') as f:
            data = f.read()
            data = data.replace(possible_old_ff_code, ligand_ff_code)
        with open(receptor_dict[receptor]["topol"], 'w') as f:
            f.write(data)

        for ligand in ligand_dict:

            # Make ligand topologies
            assemble_path = os.path.join(out_dir, receptor, ligand)
            assemble_path = os.path.abspath(assemble_path)
            ASSEMBLE_PDB = make_ligand_topology(
                receptor_code = receptor,
                ligand_code = ligand,
                docking_dict = docking_dict,
                ligand_dict = ligand_dict,
                ligand_ff_code = ligand_ff_code,
                openff_code = openff_code,
                out_dir = out_dir,)
            ASSEMBLE_PDB.atoms = pdb.PDB(receptor_dict[receptor]["update"]).atoms + pdb.PDB(receptor_dict[receptor]["membrane"]).atoms + ASSEMBLE_PDB.atoms
            # BEcasue the path is the correct in make_ligand_topology,
            # there is not need to specify it again: out_dir/receptor/ligand
            ASSEMBLE_PDB.write(os.path.join(assemble_path, "ASSEMBLY.pdb"))

            [tools.cp(item, assemble_path ,r = True) for item in [receptor_dict[receptor]["toppar"], receptor_dict[receptor]["topol"], os.path.join('CHARMM-GUI', '*.mdp')]]

            tools.cp(os.path.join('Slipids_2020.ff', "Slipids_2020_itp_files", "POPC.itp"), os.path.join(assemble_path, 'toppar'))

            #Now solvate the system
            gmxsolvate(assemble_path, "ASSEMBLY.pdb", "topol.top", receptor_dict[receptor]["cryst1"]["vector"], editconf_angles = receptor_dict[receptor]["cryst1"]["angles"], editconf_bt = "tric", solvate_cs = "spc216", out_file = "ASSEMBLY.pdb")

            #Generate the ions
            gmxions(assemble_path, "ASSEMBLY.pdb", "topol.top", out_file = "ASSEMBLY.pdb")

            #Now the index file is generated
            finalndx(assemble_path, "ASSEMBLY.pdb", chainsID = PROT_PDB.get_chains())

            GROMACS_JOB(os.path.join(assemble_path, "job.sh"), hostname=hostname, GROMACS_version=GROMACS_version)

    os.chdir(cwd)
    tmp_wd.cleanup()

if __name__ == "__main__":
    pass



