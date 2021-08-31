"""File containing functions to convert .csv format to .vtx format for OpenCARP simulations"""

import pandas as pd
import argparse
import numpy as np
from tqdm import tqdm
import glob


def convert(filename="stim", meshdir='../'):  # converts paraview csv to something openCarp can work with
    """
    Takes a csv file and writes a vtx file in the same directory
    :param filename: the name of the paraview csv file
    :return: 0 if successful
    """
    # input
    df = pd.read_csv(filename, dtype=int)

    # output
    name = filename.split(".")[0]
    outfile = open(name+".vtx", "w+")  # ID's of tagged regions
    outfile_vis = open(name+".dat", "w+")  # all pts. 1 if tagged, 0 if not. For meshalyzer.
    if "/" in filename:
        name = glob.glob("/".join(filename.split("/")[:-1]) + "/*.pts")[0]  # takes first .pts file
    else:
        name = glob.glob(meshdir+"*.pts")[0]
    with open(name) as file:
        n_points = int(file.readline())

    tagpts = np.zeros(n_points, dtype=int)  # n points of mesh, used for visualisation in meshalyzer

    # header
    outfile.write(str(len(df.index)) + '\nintra\n')  # amount and type of vertices

    for e in df.values:  # writing all ID's that are blocks
        tagpts[e[0]] = 1
        outfile.write(str(e[0])+'\n')

    for e in tagpts:  # pts_t file, where index = 1 if point belongs to paraview selected regions
        #  for visualisation in meshalyzer
        outfile_vis.write(str(e)+'\n')

    outfile.close()
    outfile_vis.close()

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # fetching the arguments
    parser.add_argument('--infiles',
                        help="Relative paths to stim and/or block file. Separate multiple files with a comma",
                        type=str, default="")
    parser.add_argument('--meshdir',
                        help="Relative paths to stim and/or block file. Separate multiple files with a comma",
                        type=str, default="../")
    args = parser.parse_args()
    if not args.infiles:
        stims = glob.glob("stim*.csv")
        blocks = glob.glob("block*.csv")
        files = stims + blocks
    else:
        files = args.infiles
    files.sort()
    t = tqdm(files)
    for f in t:
        t.set_description(f)
        convert(f, meshdir=args.meshdir)
    print("All files converted")
