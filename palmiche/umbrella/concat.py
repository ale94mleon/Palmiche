#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, tempfile, argparse
from palmiche.utils import tools
import multiprocessing as mp
import tqdm



def main(input_path, conf_file_name, output_path = '.', cpu = 0):
    """Take the last frame of the production of each window, and then concatenate them.

    Args:
        input_path (path): The path were the windows directories are with name xxxxx, where x is a number
        conf_file_name (string): the name of the last frame exported when the umbrella is finished
        output_path (str, optional): The output path where the file cat_windows.xtc will be output. Defaults to '.'.
    """
    pattern = re.compile("^\d{5}$")
    conf_file_paths = sorted([os.path.join(input_path, window, conf_file_name) for window in tools.list_if_dir(input_path) if pattern.match(window)])
    temp_folder = tempfile.TemporaryDirectory(dir = '.', prefix='.')
    #with tempfile.TemporaryDirectory() as temp_folder:
    xtcs = []
    cmds = []
    for (i, conf) in enumerate(conf_file_paths):
        xtc_tmp = os.path.join(temp_folder.name, f"{i}.xtc")

        cmds.append(f"export GMX_MAXBACKUP=-1;  gmx trjconv -f {conf} -o {xtc_tmp}")
        # I catch the xtc names and then they will be ordered
        xtcs.append(xtc_tmp)
   #In order to parallelize the dist calculation
    if cpu:
        jobs = cpu
    else:
        jobs = mp.cpu_count()
    pool = mp.Pool(jobs)
    print("Converting the last frame of each window to .xtc...")
    for i in tqdm.tqdm(pool.imap_unordered(tools.run, cmds), total=len(cmds)):
        pass
    #pool.map(tools.run,cmds)
    pool.close()
    print("Done!")

    print("Concatening .xtc files...")
    tools.run(f"gmx trjcat -f {' '.join([xtc for xtc in xtcs])} -o {os.path.join(output_path, 'cat_windows.xtc')} -cat")
    print("Finished!")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--input_path','-i',
                        help = "The path were the windows directories are with name xxxxx, where x is a number",
                        dest='input_path',
                        type=str)
    parser.add_argument('--conf_file_name', '-c',
                        help = "the name of the last frame exported when the umbrella is finished",
                        dest='conf_file_name',
                        type=str)
    parser.add_argument('--output_path', '-o',
                        help = "The output path where the file cat_windows.xtc will be output. Defaults to '.'",
                        default='.',
                        type=str)

    args = parser.parse_args()

    main(args.input_path, args.conf_file_name, args.output_path)

