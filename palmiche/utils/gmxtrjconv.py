#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from palmiche.utils import tools
import argparse
import tempfile, os

def gmxtrjconv(tpr, conf, xtc, trans = [0,0,0], dt = 0, vmd = False):
    """This function try to create a nice visualization of th MD simulation.
    The idea is test several translation till the desire part of the system is as far as possible from the edges. In this way will not break down
    through the periodic boundaries. THis "breaking" would happen for system with more that one molecule (a pentamer) when one of them goes to
    "other side" of the periodic boundary. But for a molecule itself i will not break.

    Args:
        tpr (str path-like): tpr file
        conf (str path-like): configuration file, pdb or gro extensions
        xtc (str path-like): trajectory file
        trans (list, optional): the x, y, z translation amounts. Defaults to [0,0,0].
        dt (int, optional): The interval of time to be considered. For initial testing is wise to use a high value, that give only a few frame of the trajectory. Defaults to 0.
        vmd (bool, optional): This will try to open vmd in order that you visualize the system on the testing phase. Defaults to True.
    """
    if dt:
        cmd_add = f'-dt {dt}'
    else:
        cmd_add = ''
    tmpdir = tempfile.TemporaryDirectory(dir = '.', prefix='.')


    tools.run(f"""
        export GMX_MAXBACKUP=-1
        echo "system" | gmx trjconv -s {tpr} -f {conf} -o {os.path.join(tmpdir.name, 'tmp.gro')} -trans {" ".join([str(i) for i in trans])}
        echo "system" | gmx trjconv -s {tpr} -f {os.path.join(tmpdir.name, 'tmp.gro')} -o processed.gro -pbc mol -ur compact

        echo "system" | gmx trjconv -s {tpr} -f {xtc} -o {os.path.join(tmpdir.name, 'tmp1.xtc')} -trans {" ".join([str(i) for i in trans])} {cmd_add}
        echo "system" | gmx trjconv -s {tpr} -f {os.path.join(tmpdir.name, 'tmp1.xtc')} -o {os.path.join(tmpdir.name, 'tmp2.xtc')} -pbc mol -ur compact
        echo "protein system" | gmx trjconv -s processed.gro -f {os.path.join(tmpdir.name, 'tmp2.xtc')} -fit rot+trans -o processed.xtc
        """)
    if vmd:
        tools.run("vmd processed.gro processed.xtc")



def main():
    """This is just the main function of the script, see __doc__ for details"""

    parser = argparse.ArgumentParser(description=__doc__,
                                        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s',
                        help = "The tpr file",
                        dest='tpr',
                        type=str)
    parser.add_argument('-x',
                        help = "The trajectory file, xtc",
                        dest='xtc',
                        type=str)
    parser.add_argument('-c',
                        help ='The configuration file: gro or pdb',
                        dest = 'conf',
                        type=str)
    parser.add_argument('-trans',
                        help ='The translation vector, this is how much you want ot translate the system. Specified the three numbers, defaults is 0 0 0',
                        dest = 'trans',
                        nargs = 3,
                        default = [0,0,0],
                        type=float)
    parser.add_argument('-dt',
                        help ='The interval of time to be considered. For initial testing is wise to use a high value, that give only a few frame of the trajectory. Defaults to 0.',
                        dest = 'dt',
                        default = 0,
                        type=float)
    parser.add_argument('-vmd',
                        help ='This will try to open vmd in order that you visualize the system on the testing phase. Defaults to False',
                        nargs = "?",
                        dest = 'vmd',
                        const = True,
                        default = False,
                        type=bool)
    args = parser.parse_args()



    gmxtrjconv(args.tpr, args.conf, args.xtc, trans = args.trans, dt = args.dt, vmd = args.vmd)