#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tempfile, os, tqdm, pickle, json
from openbabel import openbabel as ob
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import SaltRemover
from rdkit.Chem import rdmolops
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
import pubchempy as pcp
import multiprocessing as mp
from datetime import datetime

 
# Consider the use of io.StringIO, in cases that 


def update_compounds(references, Threshold = 80, compounds_file = 'compounds.pickle', checkpoint_file = "cids_checkpoint.pickle"):
    # data = {"Last_update": now, "References": references, "Compounds": new_cids}
    now = datetime.now()#.strftime("%Y-%m-%d %H:%M:%S")

    print(f'Similarity search. Getting CIDs:...')
    uniques = []
    for ref in tqdm.tqdm(references):
        uniques += pcp.get_cids(ref[0] ,ref[1], searchtype='similarity', Threshold = Threshold)
    uniques = set(uniques)

    # Searching for non download entrance and adding to the queue
    if os.path.isfile(checkpoint_file):
        with open(checkpoint_file, 'rb') as f:
            uniques.update(pickle.load(f)['CID'])

    # Searching for the CID already on the database
    if os.path.isfile(compounds_file):
        with open(compounds_file, 'rb') as f:
            uniques -= set(c.cid for c in pickle.load(f)['Compounds']) # Only those CID not in the data base

    # Here I could restrict to only some part of that object if the data is to big

    print(f'Updating compounds data base {compounds_file} ...')
    new_compounds = []
    non_downloaded_cids = []
    for cid in tqdm.tqdm(list(uniques)):
        # The try except block is because the possible reject of the server
        # And in that case we need to save the data.
        try:
            new_compounds.append(pcp.Compound.from_cid(cid))
        except:
            non_downloaded_cids.append(cid)
            next
    
    
    print(f"The data base {compounds_file} was successfully updated!")
    # 
    if non_downloaded_cids:
        CIDs = {"Last_update": now, "CID": non_downloaded_cids}
        with open(checkpoint_file, 'wb') as f:
            pickle.dump(CIDs, f)
    else:# All CIDs were downloaded
        if os.path.isfile(checkpoint_file):
            os.remove(checkpoint_file)
            print(f"The {checkpoint_file} is not longer needed, therefore was removed")

    # Opening the compounds in case that exists

    if os.path.isfile(compounds_file):
        with open(compounds_file, 'rb') as f:
            data = pickle.load(f)
            # Updating
            
            # Datetime
            data["Last_update"] = now
            # References
            for ref in references:
                if ref not in data["References"]:
                    data["References"].append(ref)
            # Compounds
            data["Compounds"] += new_compounds

    else:
        # Initialization of the data
        data = {"Last_update": now, "References": references, "Compounds": new_compounds}
    
    # Writting
    with open(compounds_file, 'wb') as f:          
        pickle.dump(data, f)


    return new_compounds


def obconvert(inpath, outpath):
    """Convert  molecule ussing openbabel

    Args:
        input (str, path): input molecule.
        output (str, path): must have the extention of the molecule.
    """
    in_ext = os.path.basename(inpath).split('.')[-1]
    out_ext = os.path.basename(outpath).split('.')[-1]

    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats(in_ext, out_ext)
    mol = ob.OBMol()
    obConversion.ReadFile(mol, inpath)   # Open Babel will uncompressed automatically
    #mol.AddHydrogens()
    obConversion.WriteFile(mol, outpath)

def confgen(smiles, outpath = "mol.pdbqt"):
    """Create a 3D model from smile to

    Args:
        smiles (str): a valid smiles code.
        outpath (str, optional): The output molecule with its corresponding extension. Defaults to "mol.pdbqt".
    """
    moltemp = tempfile.NamedTemporaryFile(suffix='.mol')
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    Chem.MolToMolFile(mol, moltemp.name)
    obconvert(moltemp.name, outpath)
    with open(outpath, 'r') as o:
        string = o.read()
    return string



# def lipinski_filter_pcp(compound_obj, maxviolation = 2):
#     cont = 0
#     if compound_obj.h_bond_donor_count > 5: cont += 1
#     if compound_obj.h_bond_acceptor_count > 10: cont += 1
#     if compound_obj.exact_mass > 500: cont += 1
#     if compound_obj.xlogp > 5: cont += 1
#     #if compound_obj.rotatable_bond_count > 10: cont += 1
#     if cont >= maxviolation: return False
#     return True

# def elements_filter_pcp(compound_obj):
#     vina_elements = ['Br']
#     for atom in compound_obj.atoms:
#         if atom.element not in vina_elements:
#             return False
#     return True

# def elements_filter_pcp(mol):
#     vina_elements = ['Br']
#     for atom in mol.GetAtoms():
#         if atom.GetSymbol() not in vina_elements:
#             return False
#     return True

def lipinski_filter(mol, maxviolation = 2):
    filter = {
        'HBD': {'method':Descriptors.rdMolDescriptors.CalcNumLipinskiHBD,'cutoff':5},
        'HBA': {'method':Descriptors.rdMolDescriptors.CalcNumLipinskiHBA,'cutoff':10}, 
        'wt': {'method':Descriptors.MolWt,'cutoff':500},
        'MLogP': {'method':Descriptors.MolLogP,'cutoff':5}
        #'rotbond': {'method': Descriptors.NumRotatableBonds(), 'cutoff':10}
    }
    cont = 0
    for property in filter:
        
        if filter[property]['method'](mol) > filter[property]['cutoff']:
            cont += 1
        if cont >= maxviolation:
            return False
    return True

def remove_salt(mol):
    remover = SaltRemover.SaltRemover()
    return remover.StripMol(mol, dontRemoveEverything=True)

def neutralize_atoms(mol):
    """Neutralizing Molecules¶
    Author: Noel O’Boyle (Vincent Scalfani adapted code for RDKit)  
    Source: https://baoilleach.blogspot.com/2019/12/no-charge-simple-approach-to.html
    Index ID#: RDKitCB_33
    Summary: Neutralize charged molecules by atom.
    https://www.rdkit.org/docs/Cookbook.html?highlight=neutralize

    Args:
        mol (rdkit object): A rdkit molecules

    Returns:
        rdkit molecules: the molecule after neutralization
    """
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

def largest_mol(mol):
    """
    Author: Andrew Dalke and Susan Leung
    Source: https://sourceforge.net/p/rdkit/mailman/message/36355644/ and https://github.com/susanhleung/rdkit/blob/dev/GSOC2018_MolVS_Integration/rdkit/Chem/MolStandardize/tutorial/MolStandardize.ipynb
    Index ID#: RDKitCB_31
    Summary: Select largest fragment from a molecule
    """
    mol_frags = rdmolops.GetMolFrags(mol, asMols = True)
    if len(mol_frags) == 1:
        return mol
    return max(mol_frags, default=mol, key=lambda m: m.GetNumAtoms())
    
    #largest_Fragment = rdMolStandardize.LargestFragmentChooser()
    #return largest_Fragment.choose(mol)


def compounds_filter(pcp_compound):
    pdbqt_tmp = tempfile.NamedTemporaryFile(suffix='.pdbqt')

    smiles = pcp_compound.isomeric_smiles
    compound = {'CID': pcp_compound.cid, 'isomeric_smiles_pcp': smiles}
    mol = Chem.MolFromSmiles(smiles)
    mol = remove_salt(mol)
    
    #un = rdMolStandardize.Uncharger()
    #un.uncharge(mol)
    mol = neutralize_atoms(mol)
    mol = largest_mol(mol)
    compound['rdkit_smiles'] = Chem.MolToSmiles(mol)
    if lipinski_filter(mol):
        pdbqt_string = confgen(compound['rdkit_smiles'], outpath = pdbqt_tmp.name)
        compound['pdbqt_string'] = pdbqt_string
        return compound
    else:
        return None

def receptor(pdbqt_file, box_file, receptor_name):
    with open(pdbqt_file, 'r') as f:
        pdbqt_string = f.read()
    with open(box_file, 'r') as f:
        grid_opt = json.load(f)[receptor_name]
    return {'name': receptor_name, 'pdbqt_string': pdbqt_string, 'grid_opt': grid_opt}

def ligands(compounds, cpus = 0):

    if cpus:
        jobs = cpus
    else:
        jobs = mp.cpu_count()
    pool = mp.Pool(jobs)

    data =  list(filter(None, pool.imap_unordered(compounds_filter, compounds)))
    pool.close()
    
    dictionary = dict()
    for lig in data:
        dictionary[lig['CID']] = lig
    return dictionary

def vinaconfig(receptor_path,
                ligand_path,
                grid_center,
                size,
                num_modes = 1,
                exhaustiveness = 8,
                cpu = 3):

    name = os.path.basename(ligand_path).split(".")[0]
   
    string =f"{'receptor':25} = {receptor_path}\n"\
            f"{'ligand':25} = {ligand_path}\n\n"\
            f"{'center_x':25} = {grid_center[0]}\n"\
            f"{'center_y':25} = {grid_center[1]}\n"\
            f"{'center_z':25} = {grid_center[2]}\n\n"\
            f"{'size_x':25} = {size[0]}\n"\
            f"{'size_y':25} = {size[1]}\n"\
            f"{'size_z':25} = {size[2]}\n\n"\
            f"{'num_modes':25} = {num_modes}\n"\
            f"{'exhaustiveness':25} = {exhaustiveness}\n"\
            f"{'cpu':25} = {cpu}\n\n"\
            f"{'out'} =  {name}_out.pdbqt\n"
    return string
            
def moltosvg(mol, molSize = (300,300), kekulize = True):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg #le quite lo de replace svg: por '', pq sino no funcionaba


    
    
    

# references = [
#     ('OC(CC(O)C(F)(F)C(F)(F)F)C1CCCNC1','smiles', 'cis1_bh267_m'),
#     ('COC1CC[C@@H](C(O)CC(O)C(F)(F)C(F)(F)F)C(O)C1','smiles', 'sys_MMV007839')
# ]
# AQdicionar update info al opbjecto, la fecha, las estructuras que se emplearon, el threhold, etc..

# # Checking for necessary updates
# new_cids = update_cids(references)
# new_compounds = update_compounds(new_cids)
# for compound in new_compounds:
#     # here I need to do a lot of stuffs:
#     # Sanitze the SMILES, removes sals and so,
#     # Generate a confomation
#     # Generate a pdbqt file
#     # generate the configuation file for vina
#     # run vina
#     # save the data inside the object

#     setattr(compound, 'docking', 'a dict with final SMILES used, the config file, the exhautivines used, the initial configuration files and the BE confomation, also the energy')


"""
Puedo poner un CID que tenga los CID aceptados y los rechazados, en fin, los cid que se
han anmalizado hasta la fecha, pero se le puede a;adir ese plus,
Esto te ayuda a que si los CID que estas proponiendo ya se estudiaron ni molestarse
Tengoq u obtener los SMILES de las dos moleculas que estnaos estudiando para que sean las referencias
la cosa es que tengo que hacer una funcion que me actualice los smiles,
segun las moleculas que se adicionene
la cosa es que voy a tener dos archivos, uno
.smi
que va a tener todas los smiles procesados, y como nombre el CID
Y otro archivo .sdf con los conformeros
FInalmente una carpeta con los .pdbqt que van a tener por nombre el CID
el update va a buscar el CID en smile, tambien puede buscar por smile,
En caso de no estar, delo a;adira a los res niveles.

HE ARPENDIDO UNA PILA DE TALLAS

1- Buscar por similitud en pubchem con pubchempy
lo que hacemos es descargar solo los CID porque sino al descargar tanta informacion el servidos nos va a rechazar de
por tanto:
a = pcp.get_cids('COC1=CC2=C(C=C1)C(=O)CC(O2)(C(C(F)(F)F)(F)F)O','smiles', searchtype='similarity', Threshold = 80)
a = pcp.get_cids(2917555 ,'cid', searchtype='similarity', Threshold = 80)
todos estos metodos funcionas que:
especificas el tipo a buscar ya sea smile, cid aid o lo que quieras,
pero tinees que poner el tipo que buscaste, por
luego especificas searchtype
https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest$_Toc494865580
En esa pagina explican como se manja estos parametros
el threshold es para el tanimoto
ahora bien: generamos un set de varios cids
y luego llamamos a 


"""
if __name__ == '__main__':
    pass
    v = vinaconfig('./k', 'l', (5,8,9), (50,50,50))
    print(v)
