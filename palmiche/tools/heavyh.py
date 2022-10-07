#!/usr/bin/env python
'''
Author(s): Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>, 
\nDescription: This script modifies gromacs topology files to simulate with heavy hydrogens. This
allows to decrease angular and out-of-plane motions involving hydrogens and therefore allowing to
increase the time step. To conserve the total mass, the heavy atoms connected to hydrogens will
also have a proportionally reduced mass. Example for methyl group with factor of 3 (default, see
argument descriptions at the end of the help) :\n
           H (1.008)                     H (3.024)
           |                             |
(1.008) H--C--H (1.008)  ==>  (3.024) H--C--H (3.024)
        (12.01)                       (5.9620)

(FLOAT): Atom mass\n

Another approach is to use virtual sites, however it is not possible to do the update on the GPUs
with virtual sites. Therefore, using heavy hydrogens allow to double the time step AND to use
gmx mdrun -update gpu and gain further significant speed-up.
It is atm only possible to use a factor of 4 with pdb2gmx -heavyh thus the need of this script.

If you are simulating proteins, you can try to constrain all bonds: "constraints = all-bonds" in
the .mdp file. But if the system is too big you won't be able to do the updates on the gpus.\n

Some script launch examples:\n
1. heavyh.py : a .top file is looked in current dir and all .itp files in current directory
and included in the .top will be processed as well. Files in the *posre*.itp form and itp files
containing "/" will be ignored. This has been done like this to avoid processing files like
"amber99sb-ildn.ff/forcefield.itp", files in the root force field directory and more generally to
have better controlled of what is processed where. All hydrogen masses are multiply by default
factor (i.e. 3).
2. heavyh.py -p topol.top : same as 1. but -p needed if several .top in directory.
3. heavyh.py -p topol1.itp topol2.itp: only processes the .itp files given, does not look
for .top.
4. heavyh.py -p topol.top -factor 2.5 : same as 2. but all hydrogen masses are multiply by a
factor of 2.5.
5. heavyh.py -p topol_heavyH_3.top -revert 3.024 -factor 3 -o _revertHMR : revert HMR on hydrogens
having a mass of 3.024amu (have been originally increased by a factor 3)

Ref:
1. K Anton Feenstra, Berk Hess and Herman J C Berendsen. Improving Efficiency of Large Time-Scale
Molecular Dynamics Simulations of Hydrogen-Rich Systems. Journal of Computational Chemistry, 1999.
2. Chad W. Hopkins, Scott Le Grand, Ross C. Walker and Adrian E. Roitberg. Long-time-step molecular
dynamics through hydrogen mass repartitioning. Journal of Chemical Theory and Computation, 2015.
3. Curtis Balusek et al. Accelerating Membrane Simulations with Hydrogen Mass Repartitioning.
Journal of Chemical Theory and Computation, 2019.

TODO/Notes:\n1. It is not possible to have several Hydrogens with different masses in the same
topology file yet. If those topologies are in different .itp files included in the .top, this will
be all fine though. However, for reverting back HMR with -revert flag, make sure all your hydrogens
have the same mass as the one you've given to the flag. If this is not the case, process your
topologies separately.
2. If included .itp files in .top are not in the same directory as the .top, the script will not
process it.
3. This scripts assumes that [ atoms ] block is followed by the [ bonds ] block 
(not [ pairs ], [ dihedrals ], etc.), this has been the case in all the case systems so far, but if
this leads to issues, it should be "not too hard" to change this in the script.

'''

import argparse
from glob import glob
import re
from os import path
from shutil import move
import itertools


class NoHydroError(Exception):
    """Customized exception"""
class SeveralHydroMass(Exception):
    """Customized exception"""

def file_checker(top, itp, suff):
    """This function checks for .itp (not containing 'posre' or '/') included in .top and try to
    look for them in the current directory or in the directory name of the .top file given by -p.\n
    Check if .itp/.top files have an [ atoms ] block to modifiy, else remove it from input files.\n
    Finally, builds output names for the topologies"""


    if top != []:
        print(f'Working on following .top file: {top}\n')
        regxitp = re.compile(r'(?<=#include ")(?!.*(\b\/\b|posre)).*(?=")')
        with open(top[0], 'r') as intop:
            for line in intop:
                if regxitp.search(line):
                    if regxitp.search(line).group() not in itp:
                        itp.append(regxitp.search(line).group())

        #Check if those .itp are in current directory or in the dirname of .top given by -p
        for i in itp:
            if not path.isfile(i):
                dirname_of_top = path.dirname(top[0])
                print(f"{i} included in .top not found in './', "
                      f"let's look in dirname of the -p argument (i.e. {dirname_of_top})")
                guess_itp_path = f'{dirname_of_top}/{i}'
                if path.isfile(guess_itp_path):
                    print(f'{i} found in {dirname_of_top}')
                    itp[itp.index(i)] = guess_itp_path
                else:
                    print(f"{i} included in .top also not present in {dirname_of_top}:"
                          f" exiting script")
                    raise SystemExit()

    topols = top + itp

    regx_atmblock = re.compile(r'^\[ atoms \]\n')
    topol2rm = []

    for i in topols:
        keep = 0
        with open(i, 'r') as topfile:
            for line in topfile:
                if regx_atmblock.search(line):
                    keep = 1
        if keep == 0:
            print(f'No [ atoms ] block in {i}, this file will be ignored')
            topol2rm.append(i)

    topols = [topfile for topfile in topols if topfile not in topol2rm]

    print(f'\nThe following .itp and/or .top will be used for generating heavy H topologies:'
          f'\n{topols}\n')

    topolout = [re.sub(r'(\.itp|\.top)', '', path.basename(name))
                +suff
                +'.'
                +name.split('.')[-1] for name in topols]

    topolin_topolout = {topols[i]: topolout[i] for i in range(len(topols))}

    return topolin_topolout

def top_parser(topol, revert):
    """This function is parsing the topol file and creates 2 main dictionaries: toplines and
    atomblock. The file is read line by line such that everything that is not [ atoms ] is put
    in toplines associted to their line index. All [ atoms ] lines are put in atomblock with
    associated line index. Bondlist is the lists of bonds."""

    toplines = {}
    atomblock = {}
    hydrogen_mass = set()
    hydroid2line = {}
    heavyid2line = {}
    bondlist = []
    endblock = re.compile(r"^\[ *")

    #Sanity check variable
    massin = 0

    with open(topol, 'r') as top:
        for index, atom_line in enumerate(itertools.takewhile(lambda x: '[ atoms ]' not in x, top)):
            toplines[index] = atom_line

        toplines[max(toplines.keys())+1] = '[ atoms ]\n'

        for index, atom_line in enumerate(itertools.takewhile(lambda x: '[ bonds ]' not in x, top),
                                          max(toplines.keys())+1):
            #must be before split(), because index is out of range otherwise but NO ERROR MESSAGE!
            if atom_line == '\n':
                toplines[index] = atom_line
                continue
            atom_line_list = atom_line.split()
            #print(atom_line_list)
            if atom_line_list[0] == ";":
                toplines[index] = atom_line
            elif len(atom_line_list) >= 8:
                if revert != 0:
                    if float(atom_line_list[7]) == revert:
                        hydrogen_mass.add(atom_line_list[7])
                        hydroid2line[atom_line_list[0]] = index
                        atomblock[index] = atom_line_list
                        massin += float(atom_line_list[7])
                    else:
                        heavyid2line[atom_line_list[0]] = index
                        atomblock[index] = atom_line_list
                        massin += float(atom_line_list[7])
                else:
                    if 0.99 <= float(atom_line_list[7]) <= 1.1:
                        hydrogen_mass.add(atom_line_list[7])
                        hydroid2line[atom_line_list[0]] = index
                        atomblock[index] = atom_line_list
                        massin += float(atom_line_list[7])
                    else:
                        heavyid2line[atom_line_list[0]] = index
                        atomblock[index] = atom_line_list
                        massin += float(atom_line_list[7])

        #Sanity check
        if hydrogen_mass == set():
            raise NoHydroError
        if len(hydrogen_mass) > 1:
            error = (f'Hydrogens with different masses where found in {topol}, this is not'
                     f' supported yet. Please remove outputs and improve the script.')
            raise SystemExit(error)
        #End of sanity check

        toplines[max(toplines.keys())+1] = '[ bonds ]\n'

        for index, atom_line in enumerate(top, max(toplines.keys())+1):
            if endblock.search(atom_line):
                next_block = atom_line
                break
            if atom_line == '\n':
                toplines[index] = atom_line
                continue
            atom_line_list = atom_line.split()
            if atom_line_list[0] == ";":
                toplines[index] = atom_line
            else:
                bondlist.append(atom_line.split()[0:2])
                toplines[index] = atom_line

        toplines[max(toplines.keys())+1] = next_block

        for index, atom_line in enumerate(top, max(toplines.keys())+1):
            toplines[index] = atom_line

    hydrogen_mass = list(hydrogen_mass)

    return (toplines, atomblock, float(hydrogen_mass[0]),
            hydroid2line, heavyid2line, bondlist, massin)


def atom_block_modifier(hydrogen_mass, massfactor, bondlist, hydroid2line, heavyid2line, atomblock):
    """This function iterate over the bonds columns to check which atom ids are linked to an H.
       Heavy atoms are then retrieved in atomblock thanks to heavyid2line in order to decrease
       their masses. Hydrogens atoms are also retrieved in atomblock but thanks to hydroid2line and
       their masses are increased according to the massfactor. Usually H are in 2nd col of the
       [ bonds ] block, but the elif make sure to also check col 1 just in case. The use of set
       objects makes sure to do this really fast in contrast to lists (O(n) vs O(1))."""

    new_hydrogen_mass = hydrogen_mass*massfactor
    massdelta = hydrogen_mass - new_hydrogen_mass

    #Sanity check variable
    hbond_nbr = 0

    for col1col2 in bondlist:
        if col1col2[1] in hydroid2line.keys():
            #changing heavy atm mass
            line2modify = heavyid2line[col1col2[0]]
            former_mass = float(atomblock[line2modify][7])
            new_mass = former_mass + massdelta
            atomblock[line2modify][7] = str(new_mass)
            #changing hydrogen mass
            line2modify = hydroid2line[col1col2[1]]
            atomblock[line2modify][7] = str(new_hydrogen_mass)
            hbond_nbr += 1
        elif col1col2[0] in hydroid2line.keys():
            line2modify = heavyid2line[col1col2[1]]
            #atom mass is in 8th column
            former_mass = float(atomblock[line2modify][7])
            new_mass = former_mass + massdelta
            atomblock[line2modify][7] = str(new_mass)
            #changing hydrogen mass
            line2modify = hydroid2line[col1col2[0]]
            atomblock[line2modify][7] = str(new_hydrogen_mass)
            hbond_nbr += 1

    if hbond_nbr != len(hydroid2line):
        error = ('Error: the number of parsed hydrogens in hydrogenid variable is not the same as '
                 'the number of hydrogen bonds in hbond_nbr. Something is wrong !')
        raise SystemExit(error)

    massout = sum([float(column[7]) for column in atomblock.values()])

    return atomblock, massout


def write(toplines, atomblock, out):
    """This function reformates the modified atomblock, merges the toplines and the formatted
       modified atomblock and finally writes the ouptput topologies"""

    atomblock_f = {}
    for linenbr, atmline_list in atomblock.items():
        if len(atmline_list) == 8:
            atomblock_f[linenbr] = (f'{atmline_list[0]:>7}{atmline_list[1]:>12}'
                                    f'{atmline_list[2]:>8}{atmline_list[3]:>8}'
                                    f'{atmline_list[4]:>8}{atmline_list[5]:>8}'
                                    f'{f"{float(atmline_list[6]):.11g}":>12}'
                                    f'{f"{float(atmline_list[7]):.11g}":>12}\n')
        else:
            junk = ' '+' '.join([str(elem) for elem in atmline_list[8::]])
            atomblock_f[linenbr] = (f'{atmline_list[0]:>7}{atmline_list[1]:>12}'
                                    f'{atmline_list[2]:>8}{atmline_list[3]:>8}'
                                    f'{atmline_list[4]:>8}{atmline_list[5]:>8}'
                                    f'{f"{float(atmline_list[6]):.11g}":>12}'
                                    f'{f"{float(atmline_list[7]):.11g}":>12}'
                                    f'{junk}\n')

    #merging toplines with atomblock: entire new topology file
    toplines.update(atomblock_f)

    #To do: make sanity check: all toplines.keys should be different and increasing from 0 to
    #len(toplines)
    with open(out, 'w') as outtop:
        for line_nbr in range(0, len(toplines)):
            outtop.write(toplines[line_nbr])


def run_heavyh(topol = "./*.top", massfactor = 3, revert = 0, suff = "_heavyH"):
    """
    

    Parameters
    ----------
    topol : TYPE, optional, str, path, or regular expresion
        DESCRIPTION. The default is "./*.top".
        Give input topology (.top and/or itp), default: ./*.top, 
        if several .top in\ncurrent dir or given, error will be raised.
        However you can give as many .itp as you want (+max 1 .top)
    massfactor : TYPE, optional, float
        DESCRIPTION. The default is 3. int
        Give the factor by which each hydrogen mass should be multiplied or
        divided if you use revert.
    revert : TYPE, optional, float
        DESCRIPTION. The default is 0.
        Used if you want to decrease H mass instead, you should give the H
        mass in\nyour topology as argument for this flag
    suff : TYPE, optional, str
        DESCRIPTION. The default is "_heavyH".
        Suffix for output topology file name.

    Raises
    ------
    SystemExit
        DESCRIPTION.

    Returns
    -------
    None.

    """
    topol = glob(topol)

    top = [top for top in topol if ".top" in top]
    itp = [itp for itp in topol if ".itp" in itp]

    #Manage only one .top/inputs, this simplifies the script.
    if len(top) > 1:
        error = ('Several .top files in the current directory or given, please use -p to precise'
                 ' which .top should be used (you can however give as many .itp as you want or even'
                 ' just one .itp as input if you want). Tip: if you want to process several .top,'
                 ' just use a bash for loop over your .top using this scripti.')
        raise SystemExit(error)

    if revert != 0:
        massfactor = 1/massfactor

    topolin_topolout = file_checker(top, itp, suff)

    #variable needed to NOT modifiy the name of .itp without H in the .top
    no_hydrogen = []

    for items in topolin_topolout.items():
        print(f"\nProcessing {items[0]}...\n")

        try:
            toplines, atomblock, hydrogen_mass, hydroid2line,   \
            heavyid2line, bondlist, massin = top_parser(items[0], revert)

            atomblock, massout = atom_block_modifier(hydrogen_mass, massfactor, bondlist,
                                                     hydroid2line, heavyid2line, atomblock)

            write(toplines, atomblock, items[1])

            print(f'Difference between overall mass in input '
                  f'and overall mass in output is: {massin-massout:.4f} (should be 0)'
                  f'\n(massin is: {massin:.7g}, massout is: {massout:.7g})')

        except NoHydroError:
            print(f'***CAREFUL*** No hydrogens were found in: {items[0]}. If you think it makes '
                  f'sense for this molecular entity to not have hydrogens then you'
                  f'can ignore what follows.\n'
                  f'Else possible reasons could be:\n1. The hydrogen masses in the input are not '
                  f'between 0.99 and 1.1, this means that the script needs to be optimized.'
                  f'\n2. Because the top formatting is unusual'
                  f', you should review the code in this case too.\n')
            no_hydrogen.append(items[0])

    #This replaces the name of the modified .itp in the .top file if .itp contains H
    if top != []:
        #case where .top does not contain H
        if top[0] in no_hydrogen:
            new_name = top[0].replace('.top', suff+'.top')
            with open(top[0], 'r') as intop, open(new_name, 'w') as outtop:
                for line in intop:
                    for i in topolin_topolout:
                        if i not in no_hydrogen and i in line:
                            line = line.replace(i, topolin_topolout[i])
                    outtop.write(line)
        #case where .top does not contain atom block
        if top[0] not in topolin_topolout:
            new_name = top[0].replace('.top', suff+'.top')
            with open(top[0], 'r') as intop, open(new_name, 'w') as outtop:
                for line in intop:
                    for i in topolin_topolout:
                        if i not in no_hydrogen and i in line:
                            line = line.replace(i, topolin_topolout[i])
                    outtop.write(line)
        #case where .top does contain H
        else:
            with open(topolin_topolout[top[0]], 'r') as outtop, open('temp.top', 'w') as temp:
                for line in outtop:
                    for i in topolin_topolout:
                        if i not in no_hydrogen and i in line:
                            line = line.replace(i, topolin_topolout[i])
                    temp.write(line)
            #print('moving ./temp.top to', "./"+topolin_topolout[top[0]])
            move('./temp.top', "./"+topolin_topolout[top[0]])    



# If called from the command line...
if __name__ == "__main__":
    
    """This is just the main function of the script executing top_parser, atom_block_modifier
    and write functions, see __doc__ for details"""

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p',
                        default=glob('./*.top'),
                        help='Give input topology (.top and/or itp), default: ./*.top, '
                             'if several .top in\ncurrent dir or given, error will be raised.'
                             ' However you can give as many .itp as you want (+max 1 .top)',
                        dest='topol',
                        nargs='*')
    parser.add_argument('-factor',
                        default=3,
                        help='Give the factor by which each hydrogen mass should be multiplied or'
                             ' divided if you use the -revert flag'
                             ' (default: %(default)s)',
                        dest='massfactor',
                        type=float)
    parser.add_argument('-revert',
                        default=0,
                        help='Used if you want to decrease H mass instead, you should give the H'
                             ' mass in\nyour topology as argument for this flag',
                        dest='revert',
                        type=float)
    parser.add_argument('-o',
                        default='_heavyH',
                        help=f'Suffix for output topology file name (default: %(default)s)',
                        action='store',
                        dest='suff',
                        type=str)
    args = parser.parse_args()

    top = [top for top in args.topol if ".top" in top]
    itp = [itp for itp in args.topol if ".itp" in itp]

    #Manage only one .top/inputs, this simplifies the script.
    if len(top) > 1:
        error = ('Several .top files in the current directory or given, please use -p to precise'
                 ' which .top should be used (you can however give as many .itp as you want or even'
                 ' just one .itp as input if you want). Tip: if you want to process several .top,'
                 ' just use a bash for loop over your .top using this scripti.')
        raise SystemExit(error)

    if args.revert != 0:
        args.massfactor = 1/args.massfactor

    topolin_topolout = file_checker(top, itp, args.suff)

    #variable needed to NOT modifiy the name of .itp without H in the .top
    no_hydrogen = []

    for items in topolin_topolout.items():
        print(f"\nProcessing {items[0]}...\n")

        try:
            toplines, atomblock, hydrogen_mass, hydroid2line,   \
            heavyid2line, bondlist, massin = top_parser(items[0], args.revert)

            atomblock, massout = atom_block_modifier(hydrogen_mass, args.massfactor, bondlist,
                                                     hydroid2line, heavyid2line, atomblock)

            write(toplines, atomblock, items[1])

            print(f'Difference between overall mass in input '
                  f'and overall mass in output is: {massin-massout:.4f} (should be 0)'
                  f'\n(massin is: {massin:.7g}, massout is: {massout:.7g})')

        except NoHydroError:
            print(f'***CAREFUL*** No hydrogens were found in: {items[0]}. If you think it makes '
                  f'sense for this molecular entity to not have hydrogens then you'
                  f'can ignore what follows.\n'
                  f'Else possible reasons could be:\n1. The hydrogen masses in the input are not '
                  f'between 0.99 and 1.1, this means that the script needs to be optimized.'
                  f'\n2. Because the top formatting is unusual'
                  f', you should review the code in this case too.\n')
            no_hydrogen.append(items[0])

    #This replaces the name of the modified .itp in the .top file if .itp contains H
    if top != []:
        #case where .top does not contain H
        if top[0] in no_hydrogen:
            new_name = top[0].replace('.top', args.suff+'.top')
            with open(top[0], 'r') as intop, open(new_name, 'w') as outtop:
                for line in intop:
                    for i in topolin_topolout:
                        if i not in no_hydrogen and i in line:
                            line = line.replace(i, topolin_topolout[i])
                    outtop.write(line)
        #case where .top does not contain atom block
        if top[0] not in topolin_topolout:
            new_name = top[0].replace('.top', args.suff+'.top')
            with open(top[0], 'r') as intop, open(new_name, 'w') as outtop:
                for line in intop:
                    for i in topolin_topolout:
                        if i not in no_hydrogen and i in line:
                            line = line.replace(i, topolin_topolout[i])
                    outtop.write(line)
        #case where .top does contain H
        else:
            with open(topolin_topolout[top[0]], 'r') as outtop, open('temp.top', 'w') as temp:
                for line in outtop:
                    for i in topolin_topolout:
                        if i not in no_hydrogen and i in line:
                            line = line.replace(i, topolin_topolout[i])
                    temp.write(line)
            #print('moving ./temp.top to', "./"+topolin_topolout[top[0]])
            move('./temp.top', "./"+topolin_topolout[top[0]])
