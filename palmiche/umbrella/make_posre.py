#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from palmiche.utils import tools, pdb, gro
import os

                
                    
def union(list_of_list):
    """
    

    Parameters
    ----------
    list_of_list : TYPE, list
        DESCRIPTION. like
        [[1,2,1,3,4,5],[5,5,6,6,4,1,8,1],[1,1,10,5,20,36],[55,5,1,45,1]]


    Returns
    -------
    TYPE
        a list with the elemnts in comun between all the lists

    """
    to_return = list_of_list[0]
    
    for item in list_of_list:
        to_return = set(to_return) & set(item)  
    #    to_return = list(filter(lambda x: x in to_return, item))
    return sorted(to_return)
    #return sorted(list(set(to_return)))

def ndx_rectification(index_file, config_path, grouplist):
    """
    This function is for the case that you have more than one monomer fo the 
    same protein. In this case each monomer should have the same residues (number)
    This function looks for the common residues of the groups decleared in 
    grouplist and then look for the index atom of these aminoacids 
    on each monomer. Then it will try to guess to which group belong the new
    index and finally export a dictionary that could be write to a GROMACS
    index file using the function the method write of the class NDX.
    The new atom index included all the atoms in the residues (also Hydrogens)

    Parameters
    ----------
    index_file : TYPE
        DESCRIPTION.
    config_path : TYPE
        DESCRIPTION.
    grouplist : TYPE
        DESCRIPTION.

    Returns
    -------
    new_index_dict : TYPE
        DESCRIPTION.

    """
    old_index_dict = tools.NDX(index_file).data
    new_index_dict = old_index_dict.copy()
    basename = os.path.basename(config_path)
    ext = basename.split('.')[-1]
    to_work = []
    
    if ext == 'gro':
        atoms = gro.GRO(config_path).atoms
        
        for group in grouplist:
            residues = []
            for old_index in old_index_dict[group]:
                residues.append(atoms[old_index-1].resNumber)
            to_work.append(list(set(residues)))

        common_resi = union([i for i in to_work])
        print("These are the aa taken into account, use the following command to see them in Pymol:\n show sticks, resi "+"+".join([str(i) for i in common_resi]))
        
        
        test = atoms[0].resNumber
        tmp_index = []
        total_index = []
        for i, atom in enumerate(atoms):
            if atom.resNumber < test or i == len(atoms)-1:
                    
                total_index.append(tmp_index)
                tmp_index = []                
            if atom.resNumber in common_resi:
                tmp_index.append(i+1)
            test = atom.resNumber

        #Guesing the groups
        for group in grouplist:
            old_index = new_index_dict[group]
            intersection = 0

            for new_index in total_index:
                tmp = len(set(old_index) & set(new_index))
                
                if tmp > intersection:
                    intersection = tmp
                    new_index_dict[group] = new_index
            print(f"The new index group of {group} has {intersection} atoms that also were in the old index group that had a total of {len(old_index)} atoms")
        
        return new_index_dict
    
    elif ext == 'pdb':
        atoms = pdb.PDB(config_path).atoms
        
        for group in grouplist:
            residues = []
            for old_index in old_index_dict[group]:
                residues.append(atoms[old_index-1].resSeq)
            to_work.append(list(set(residues)))

        common_resi = union([i for i in to_work])
        print("These are the aa taken into account, use the following command to see them in Pymol:\n show sticks, resi "+"+".join([str(i) for i in common_resi]))
        
        
        test = atoms[0].resSeq
        tmp_index = []
        total_index = []
        for i, atom in enumerate(atoms):
            if atom.resSeq < test or i == len(atoms)-1:
                    
                total_index.append(tmp_index)
                tmp_index = []                
            if atom.resSeq in common_resi:
                tmp_index.append(i+1)
            test = atom.resSeq

        #Guesing the groups
        for group in grouplist:
            old_index = new_index_dict[group]
            intersection = 0

            for new_index in total_index:
                tmp = len(set(old_index) & set(new_index))
                
                if tmp > intersection:
                    intersection = tmp
                    new_index_dict[group] = new_index
            print(f"The new index group of {group} has {intersection} atoms that also where in the old index group that had a total of {len(old_index)} atoms")
        
        return new_index_dict
    
    else:
        return(f"The file {config_path} has not a gro or pdb extension")

                  

def flat_bottom_itp_maker(file_path, geometry = 8, r = 'FLAT_BOTOM_r', k = 'FLAT_BOTOM_k', output = 'flat_bottom_posre.itp', H_atoms = False, backup = True):
    """
    

    Parameters
    ----------
    file_path : TYPE path to a gro pdb or itp file
        The configuration fiel from where .
    geometry : TYPE (int), optional
        DESCRIPTION. The default is 8.
    r : TYPE (float)
        DESCRIPTION. The default is 'FLAT_BOTOM_k'. And this flat is controlled with the mdp. You could set a number also, 1.0
    output : TYPE (string), optional
    k : TYPE (str or float), optional
        DESCRIPTION. The default is 'FLAT_BOTOM_k'. And this flat is controlled with the mdp. You could set a number also, 600
    output : TYPE (string), optional
        DESCRIPTION. The default is 'flat_bottom_posre.itp'.
    H_atoms : TYPE (bool), optional
        DESCRIPTION. The default is False. Control the flag of the function
        tools.get_atom_index()
    
    For the description of geometry, r and k, please consult the GROMACS manual:
        https://manual.gromacs.org/
    Returns
    -------
    The posre file with the flat bottom potential option
    if and itp file was given, it will add at the end of the itp the
    corresponding reference lines for the posre file.

    """
    func_type = 2 #This is for flat bottom potential in GROMACS
    if backup: tools.backoff(output)
    atom_index = tools.get_atom_index(file_path, H_atoms=H_atoms)
    
    basename = os.path.basename(file_path)
    ext = basename.split('.')[-1].lower()
    to_print = [f" {index:<15}{func_type:<15}{geometry:<15}{r:<15}{k:<15}\n" for index in atom_index]
    with open(output,'w') as out:
        out.write(";Generated by make_posre.py\n"\
                 f";{'index':<15}{'func_type':<15}{'geometry':<15}{'r':<15}{'k':<15}\n")
        for item in to_print:
            out.write(item)
    if ext == 'itp':
        if backup: tools.backoff(file_path)
        with open(file_path, 'a') as f:
            f.write("\n; Include Position restraint file for Flat Bottom Potential\n"\
                    "#ifdef POSRES\n"\
                    f"#include \"{output}\"\n"\
                    "#endif\n")
    
def flat_bottom_config_posre(config_path,resNameList, H_atoms= False, backup = True):
    """
    !!!!!Tengo que mejorar esta funcion, creo que es mejor trabajr con el index file
    tanto para este como para  el itp generator, porque de este modo, se puede ser mas general
    y la seleccion no se restringe a solo una molecula
    Ademas, esto dara problemas en el caso que se le aplique a una seccion de la proteina
    y como mi caso, sea un penimport gro_reader, pdb_readertamero,
    pues aplicara el potencial en todos los monomeros,
    debido a que los aa se llaman iguales,
    y tal bez no es lo que se desee!!!!!!

    Parameters
    ----------
    config : TYPE path 
        DESCRIPTION. configuration file (pdb or gro)
    resNameList : TYPE list
        DESCRIPTION. List of Residues names for which the flatt
        bottom will be used
        ['LIA', 'LIB', 'LIC', ...]
    H_atoms : TYPE bool
        DESCRIPTION. The default is False. Control the flag of the function
            mol_selection_center(), 

        

    Returns
    -------
    The function will take the center for each residues in resNameList
    And then the xyz values of                 else:
                    select.append(atom)the instance will be changes for the selection center
    After this, the config file will be updated, and the original will be backup

    """
    basename = os.path.basename(config_path)
    ext = basename.split('.')[-1]
    if ext == 'gro':
        config_data = gro.GRO(config_path)
        for resName in resNameList:
            xc, yc, zc = config_data.centroid(resName_selection = [resName], H_atoms=H_atoms)
            for atom in config_data.atoms:
                if atom.resName == resName:
                    if H_atoms:
                        atom.x = xc
                        atom.y = yc
                        atom.z = zc
                    else:
                        if atom.atomName[0] not in 'Hh':
                            atom.x = xc
                            atom.y = yc
                            atom.z = zc                        
                
    elif ext == 'pdb':
        config_data = pdb.PDB(config_path)
        for resName in resNameList:
            xc, yc, zc = config_data.centroid(resName_selection = [resName], H_atoms=H_atoms)
            for atom in config_data.atoms:
                if atom.resName == resName:
                    if H_atoms:
                        atom.x = xc
                        atom.y = yc
                        atom.z = zc                       
                    else:
                        if atom.name[0] not in 'Hh':
                            atom.x = xc
                            atom.y = yc
                            atom.z = zc
    else:
        return('The extension file is not pdb or gro')
    config_data.write(backup = backup)
        


if __name__ == '__main__':
    pass
