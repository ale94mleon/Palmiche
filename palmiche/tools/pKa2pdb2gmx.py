#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
Created on    : Fri Feb 26 10:10:40 2021
Author        : Alejandro Martínez León
Mail          : [alejandro.martinezleon@uni-saarland.de, ale94mleon@gmail.com]
Affiliation   : Jochen Hub's Biophysics Group
Affiliation   : Faculty of NS, University of Saarland, Saarbrücken, Germany
===============================================================================
DESCRIPTION   :
DEPENDENCIES  :
===============================================================================
"""

from palmiche.tools import tools, pdb
import os, tempfile, argparse




'''
For pdb2gmx by default:
     LYS is protonated
     HIS, ASP and GLU deprotonated.
In case of:
    pKa > pH [HA] > [A-], it is protonated
    pKa < pH [HA] < [A-], it is deprotonated
The script will check this from the propka calculation and if the default option 
need to be changed, this will do automatically executing pdb2gmx with the command at the end:
This use propka3 not propka as command line.
You can specifie how you call propka on your computer with the flag --pka\n


From pdb2gmx

Which LYSINE type do you want for residue 14
0. Not protonated (charge 0) (LYN)
1. Protonated (charge +1) (LYS)

Which GLUTAMIC ACID type do you want for residue 3
0. Not protonated (charge -1) (GLU)
1. Protonated (charge 0) (GLH)

Which ASPARTIC ACID type do you want for residue 84
0. Not protonated (charge -1) (ASP)
1. Protonated (charge 0) (ASH)

Which HISTIDINE type do you want for residue 18
0. H on ND1 only (HID)
1. H on NE2 only (HIE)
2. H on ND1 and NE2 (HIP)
3. Coupled to Heme (HIS1)


Nomenclature of the proton in the pdb after pdb2gmx
HD1 ----> proton on delta (0 option)
HE2 ----> proton on epsilon (1 option) (the most away proton from the R chain)
could be both protons at the same time.
'''
 
def check_protonation(aa, pKa, pH):
    '''
    
    Parameters
    ----------
    aa : String
        The aa code, only ASP, GLU, HIS, LYS, case sensitive.
    pKa : float
        DESCRIPTION.
    pH : TYPE
        DESCRIPTION.

    Returns
    -------
    protonation_state : int
        0 or 1 for deprotonated or protoneted (1 and 2 for HIS).
    check : Bool
        True if the protonation match with the default of pdb2gmx.

    '''
    # 0 and 1 match with the input for pdb2gmx, except for HIS
    # The protonated state is 2 and the deprotonated state that we will use is
    # 1. H on NE2 only (HIE)
    # The next dict is only for check the protonation state. In HIS the default is deprotoneted but not precisilly 0. 
    # When deprotonated, gmx will try to figure out which state is presented (delta or epsilon)
    # But for convenience, we could use 0 and one to check if it is protoneated or deprotoneted and 
    # at the end sum 1. That means that this will not give nor the delta protonation neither Coupled to Heme (HIS1) states of HIS only
    # epsilon and protonated. This could be improved in the future and get more information from propka. For our propouse it is sufficient.
    # for the rest of the aa.
    
    gmx_default = {
        'ASP': 0,    #deprotonated
        'GLU': 0,    #deprotonated
        'HIS': 0,    #deprotonated
        'LYS': 1,    #protonated
        }
    
        
    protonation_state = 0    
    check = True
    
    if pKa > pH:
        protonation_state = 1
    if gmx_default[aa] != protonation_state:
        check = False
    if aa == 'HIS': protonation_state += 1 #(1---->H on NE2 only (HIE), 2----> dobleprotonated )
    
    return protonation_state, check



def pka2gmx(pdb_path,
            protein_forcefield = "amber99sb-ildn",
            water_forcefield = "tip3p",
            ignh = True,
            ph = 7.0,
            pKa = "propka3",
            addph = [],
            addresifile = [],
            output = "system_GMX.pdb",
            pKa_file_out = True,
            itp_top_out = True,
            num_file = True):
    # Getting cwd and changing to temporal and storaging names
    cwd = os.getcwd()
    tmpdir = tempfile.TemporaryDirectory()
    pdb_path = os.path.abspath(pdb_path)
    output = os.path.abspath(output)
    file_name = f"{os.path.basename(pdb_path).split('.')[0]}"
    pKa_file = f"{file_name}.pka"

    # Totally inefficient but needed in case of force fields in the directory
    for item in tools.list_if_dir(cwd):
        if item.endswith('.ff'):
            # I need to copy everything and not just the force field because some times you need another files outside of the ff
            # The only directories that I will copy are the ones with ends in .ff
            
            [tools.cp(file, tmpdir.name, r = True) for file in tools.list_if_file()]
            [tools.cp(dir, tmpdir.name, r = True) for dir in tools.list_if_dir() if dir.endswith('.ff')]
            break
    os.chdir(tmpdir.name)

 
    # This part if for multiple pH on the protein. Ask for the pH to consider and the files where the resi number are.

    aa_ph_manual_input = {}
    if addph and len(addph) == len(addresifile): # If you specified more than one pH and each pH has one and only resifile

        print(f'The main ph (1) is {ph}')    
        for (ph_i, resifile_i) in zip(addph,addresifile):

            with open(resifile_i, 'r') as f:
                resifile_i_lines = f.readlines()
            # Here I am assigning to each aa the corresponded specified ph. I am inverted the input, the normal input is the ph and a file with the aa that we want at that pH
            # Now I change to all the aa in the file I assign the pH, if the next file with the next pH has the same aa the pH value will be update because aa_ph_manual_input is a dictionary. 
            for resi in resifile_i_lines:
                aa_ph_manual_input[resi] = float(ph_i)

    if os.path.exists(pKa_file):
        print(f"Propka was not call because the file {pKa_file} exist.")
    else:
        tools.run(f"{pKa} {pdb_path}")


    with open(pKa_file, 'r') as f:
        pKa_lines = f.readlines()



    togmx = {
            'ASP': [],
            'GLU': [],
            'HIS': [],
            'LYS': [],
            }
    """
    In this section we are looking for the propka prediction and constructing 
    togmx: a dict of dicts that contains all the information
    {"resi_id":resi_id,"chain":chain,"prot_state":protonation_state,"check":check}
    """
    for i in range(len(pKa_lines)):
        if 'SUMMARY OF THIS PREDICTION' in pKa_lines[i]:
            j = i + 2
            while('---------') not in pKa_lines[j]:
                aa = pKa_lines[j][3:6]
                resi_id = int(pKa_lines[j][6:10])
                chain = pKa_lines[j][11]
                pKa = float(pKa_lines[j][13:21])

                if aa in togmx:
                    if resi_id in aa_ph_manual_input: # If the aa was specified with some additional pH, then check with this one.
                        protonation_state, check = check_protonation(aa, pKa, aa_ph_manual_input[resi_id])
                        print(f"The resi {'-'.join([aa, str(resi_id), chain])} was tested with the pH {aa_ph_manual_input[resi_id]}\n")
                    else:
                        protonation_state, check = check_protonation(aa, pKa, ph)
                    togmx[aa].append({"resi_id":resi_id, "chain":chain, 
                                      "prot_state":protonation_state, "check":check})

                j += 1
            break

    '''
    Propka order the aa by chain and by resi number (the same as pdb2gmx)
    
    Only we need to be sure to ordered the sequence of comands to pdb2gmx
    
    The order that pdb2gmx ask you is by alphabetic order in the case that more than one
    flag was needed. For example, if on GLU need to be fixed and two ASP
    pdb2gmx will ask you first for ASP (the chain A, all the residues of chain A
    ordered , then chain B, all residues ordered respect to the residue number). 
    Then it will ask you for GLU in the same manner
    
    Because pdb2gmx use some optimization of the H+-bond network. In order to take
    advantages of this. We check which protonation found pdb2gmx and used.
    Only if HIS give False in the check
    '''
    tools.run(f"gmx pdb2gmx -f {pdb_path} -o tmp.pdb -ff {protein_forcefield} -water {water_forcefield} -ignh -merge all")
    [tools.rm(item) for item in ["*.top","*.itp"]]
    #tools.run('grep HIS tmp.pdb | egrep "HD1|HE2" > tmp.pdb.tmp; mv tmp.pdb.tmp tmp.pdb')
    #Here I am getting all the HISs fo my system
    # This is force field depended, because the HIS could have other names in different force fields, the function fix_pdb of assemble_MD could help for it
    # This will only work for amber, but could, of course be extended.
    protonation_HIS = []
    for resi in pdb.PDB('tmp.pdb').get_residues():
        if resi.resName in ["HIS", "HIE","HIP"]: #Any possible nomenclature of HIS in the AMBER force field except for HIS1
            resatoms = [atom.name for atom in resi.atoms]
            if not set(['HD1', 'HE2']) - set(resatoms): # This means that both atoms are in resatoms
                protonation_HIS.append(2)
            elif 'HE2' in resatoms:
                protonation_HIS.append(1)
            elif 'HD1' in resatoms:
                protonation_HIS.append(0)
            else:
                pass

    if len(protonation_HIS) == len(togmx['HIS']):
        for i in range(len(togmx['HIS'])):
            # This doble check means that we only consider the non-protoneted
            #HIS to catch if is epsilon or delta. But in the case of the 
            #prediction of protonation we will keep the Propka prediction.
            #If appear a case when pdb2gmx predict protonationi and Propka did
            #not. We will keep the epsilon deprotonation state.
            if togmx['HIS'][i] != 2 and protonation_HIS[i] != 2: 
                togmx['HIS'][i]["prot_state"] = protonation_HIS[i] # only use the information if the proton is delta or epsilon
    else:
        print('''WARNING: Something terrible wrong hapened with the pdb2gmx comprobation of HIS.
              Therefore the deprotonated HIS will be considered as the epsilon only.
              ''')
    tools.rm('tmp.pdb')

    file = open('num.txt', 'wt')
          
    print('The field num.txt have the protonation state for\n')
    aa_flag_cmd = ""

    # This is for the construction of the comand and to know which aa need to be change
    aa2change = []
    ID_chains = set()
    resi_number = set()
    
    for aa in sorted(togmx): # in this way we sorted the keys and aa2change will append the aa in a sorted way
        check = set()
        for resi in togmx[aa]:
            ID_chains.add(resi["chain"])
            resi_number.add(resi["resi_id"])
            check.add(resi["check"])
            
        if False in check:
            aa_flag_cmd += f"-{aa.lower()} "
            aa2change.append(aa)
    
    ID_chains = sorted(ID_chains) #When we sorted a set, automatically you get a list
    resi_number = sorted(resi_number)
    
    # It is important the order of num.txt. for gmx (chain alphabetic--> aa alphabetic --> ID of each aa)

    for id_c in ID_chains: #(A, B, C, ...)
        for aa in aa2change: #(ASP, GLU, HIS, LYS), only for those aa that we cheked that need change before
            for number in resi_number:# (1, 2, 3, 4, ...)
                for resi in togmx[aa]: # Loops for each residue of each aa (ASP, GLU, HIS, LYS)
                    if resi["chain"] == id_c and resi["resi_id"] == number: #Here we are matching the residue with Chain id_c, aminoacid name aa, and residue number resi_id
                        print(f"{aa}-{resi}-{id_c} \t\t {resi['prot_state']}")
                        file.write(f"{resi['prot_state']}\n")
                        
    file.close() # It is important to close the file before cat it. 
    
    
    cmd = f"export GMX_MAXBACKUP=-1; cat num.txt | gmx pdb2gmx -f {pdb_path} -o {output} -ff {protein_forcefield} -water {water_forcefield} {aa_flag_cmd}"
    if ignh: cmd += " -ignh" #Just for strange cases that you are not able to remove the H of you system
    tools.run(cmd)
    print('The command used:')
    print(cmd)

    os.chdir(cwd)

    # Coping from temporal directory
    if pKa_file_out:
        tools.cp(os.path.join(tmpdir.name, pKa_file), '.')
    if itp_top_out:
        [tools.cp(item, '.') for item in [os.path.join(tmpdir.name, "*.top"), os.path.join(tmpdir.name, "*.itp")]]
    if num_file:
        tools.cp(os.path.join(tmpdir.name, 'num.txt'), '.')
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='pKa2pdb2gmx.')
    parser.add_argument("pdb_path", 
                        help = "The path to the pdb file to process.",
                        type = str)
    parser.add_argument("-pf", 
                        dest="protein_forcefield",
                        default = "amber99sb-ildn",
                        help = "The protein force field to use with GROMACS",
                        type = str)
    parser.add_argument("-wf", 
                        dest="water_forcefield",
                        default = "tip3p",
                        help = "The water force field to use with GROMACS",
                        type = str)   
    parser.add_argument("-ph", 
                        dest="ph",
                        default = 7.0,
                        help = "The general value of pH.",
                        type = float)
    parser.add_argument("-pKa", 
                        dest="pKa",
                        default = "propka3",
                        help = "The propka executable name.",
                        type = str)
    parser.add_argument("-ignh", 
                        dest="ignh",
                        nargs="?",
                        const=True,
                        default=False,
                        type=bool)
    parser.add_argument("-addph", 
                        dest="addph",
                        nargs = "+",
                        default = [],
                        help = "The other pH to condsider for each value must be file with aa to take into account with the flag addresifile.",
                        type = float)  
    parser.add_argument("-addresifile", 
                        dest="addresifile",
                        nargs = "+",
                        default = [],
                        help = "The other pH to condsider for each value must be file with aa to take into account with the flag addresifile.",
                        type = str)
    parser.add_argument('-output',
                        help='The output path, Default is system_GMX.pdb',
                        dest='output',
                        default = "system_GMX.pdb",
                        type = str) 
    parser.add_argument('-pKa_file_out',
                        help='If you want to output the pKa file. Default is False',
                        nargs="?",
                        dest='pKa_file_out',
                        const=True,
                        default=False,
                        type=bool)                   
    parser.add_argument('-itp_top_out',
                        help='If you want to output the itp and top files. Default is False',
                        nargs="?",
                        dest='itp_top_out',
                        const=True,
                        default=False,
                        type=bool) 
    parser.add_argument('-num_file',
                        help='If you want to output the file num file. This file contains the info of the protonation used for pdb2gmx. Default is False',
                        nargs="?",
                        dest='num_file',
                        const=True,
                        default=False,
                        type=bool)

    args = parser.parse_args()
    pka2gmx(args.pdb_path,
            protein_forcefield = args.protein_forcefield,
            water_forcefield = args.water_forcefield,
            ph = args.ph,
            pKa = args.pKa,
            ignh = args.ignh,
            addph = args.addph,
            addresifile = args.addresifile,
            output = args.output,
            pKa_file_out = args.pKa_file_out,
            itp_top_out = args.itp_top_out,
            num_file = args.num_file)


