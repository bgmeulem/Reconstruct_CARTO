"""File specific to this thesis. Iterates over all the carto files and tries to tetrahedralise them,
refine them or convert them.
Keeps track of which ones succeeded and which ones didn't."""

import reconstruct_mesh
import mat2csv
import os
import glob
from subprocess import Popen, PIPE
import pts2paraview
from tqdm import tqdm
import argparse
import pyvista as pv


def run(action='generateVTK', keep_intmed=False, switches="-pYkANEFq1.5/20a", HULK=False, mat=False):
    """Loops all directories in meshes/, and either
     1. creates .vtk file with default tetgen switches (generateVTK option for action)
     2. converts to carp_txt format (convert option for action)
     3. opens in paraview (paraview option for action)"""
    # switches dropped: -a2.0e+08. Tetrahedra look good and meshes run 'till the end
    if action == 'mat2ply':
        mat = True

    if mat:
        os.chdir('mat_meshes/')
    else:
        os.chdir('carto_meshes/')
        # carto meshes
        succeeded = ['OC16', 'OC19', 'OC21', 'OC22', 'OC29', 'OC32', 'OC33', 'OC34', 'OC35', 'OC36',
                     'OC37', 'OC4', 'OC41', 'OC42', 'OC45', 'OC46', 'OC49', 'OC5', 'OC50', 'OC52',
                     'OC53', 'OC56', 'OC59', 'OC63', 'OC64']
        failed = ['OC2', 'OC20', 'OC30', 'OC31', 'OC43', 'OC44', 'OC48', 'OC51', 'OC54', 'OC57',
                  'OC58', 'OC60']
        # OC2, OC31, OC44, OC51, OC57 : random floating tets, small surface of 'added points'
        # OC20: wrong hole or weird interior
        # OC43: Clusterfuck of a mesh. Looks like swiss cheese. Also super thick walls (tiny atrium?)
        # OC54: segmentation fault (interpreting string as int in reading process)
        # OC60: good mesh, but weird interior structure of heart

    dirs = os.listdir()  # all mesh directories, carto or mat
    dirs.sort()

    if action == "generateVTK":
        dirs = [e for e in dirs if e != 'Failed_Mesh_Reconstructions']

    elif not mat:
        dirs = succeeded  # only consider succeeded meshes for further actions
        tbar = tqdm(range(len(succeeded)), position=0, desc=action)

    for i, dir in enumerate(dirs):

        os.chdir(dir)  # jump in mesh directory

        if mat:  # mat format
            meshfiles = glob.glob("*.mat")
        else:  # carto format
            meshfiles = glob.glob("*.mesh")

        if action != "generateVTK" and action != "mat2ply":
            tbar.update(1)  # ignore tbar during VTK generation (output gets too messy)
        else:
            print("########################################################")
            print("#####  Directory {}/{}: {}".format(i + 1, len(dirs), dir))
            print("#####  Meshfiles: ", meshfiles)

        # Loop over all valid .mesh files in carto directories
        for j, f in enumerate(meshfiles):
            meshname = f.split(".")[0]

            # Actions
            if action == 'generateVTK':  # generate .vtk file
                print("#####  File {}/{}: {}".format(j + 1, len(meshfiles), f))

                reconstruct_mesh.run(f, switches=switches, refine_steps=10, plt=False, save=False,
                                     keep_intmed=keep_intmed, HULK=HULK)

                print("\n\t{}/{} done.\n\n".format(dir, meshname))

            elif action == 'paraview':
                # TODO: create paraview mesh, select stim, block and scars, apply and convert stimblocks to vtx
                pts2paraview.convert(meshname + ".pts")
                # subprocess.call("paraview", "{}.csv".format(meshname))

            elif action == 'convert':
                Popen("meshtool convert -imsh={0}.1 -ifmt=vtk "
                      "-omsh={0} -ofmt=carp_txt".format(meshname), shell=True, stdout=None)

            elif action == "mat2ply":
                mat2csv.mat2ply(f)

            elif action == "cleanquality":
                # threshold 0.30
                p = Popen("meshtool clean quality -msh={} -thr=.3", shell=True, stdout=PIPE, stderr=PIPE)

                # threshold 0.25
                p = Popen("meshtool clean quality -msh={} -thr=.25", shell=True, stdout=PIPE, stderr=PIPE)

                # threshold 0.20
                p = Popen("meshtool clean quality -msh={} -thr=.2", shell=True, stdout=PIPE, stderr=PIPE)

                # threshold 0.10
                p = Popen("meshtool clean quality -msh={} -thr=.1", shell=True, stdout=PIPE, stderr=PIPE)

            elif action == "writelengthtxt":
                file = open(f, 'r')
                n = 9  # line where node info is, in header
                for k, line in enumerate(file):
                    if k+1 == n:  # i starts at 0, n starts at 1
                        n_points = line.split('=')[1]
                        break
                if not os.path.exists('output'):
                    os.makedirs('output')
                with open("output/length.txt", 'w+') as of:
                    of.write(str(n_points))

        os.chdir("..")  # jump out of specific mesh directory

    os.chdir('..')  # go back to main directory, going out of meshes/ directory


if __name__ == '__main__':
    """This script loops over all folders present in directory mat_meshes/ or carto_meshes/
    Generates vtk files if action = generateVTK
    Can delete all intermediate files (so no .txt, .smesh or .csv files are saved)
    And runs meshtool convert to get all carp files needed """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    # fetching the arguments
    parser.add_argument('action',
                        type=str, default='generateVTK', choices=['generateVTK', 'convert', 'paraview', 'mat2ply',
                                                                  'cleanquality', 'writelengthtxt'],
                        help="Which action to take when looping meshes.\n"
                             "generateVTK    Tetrahedralizes each carto .mesh file. Deletes intermediate files\n"
                             "               unless --ki (Keep Intermediate) is passed. Creates .vtk\n"
                             "convert        Converts all .vtk files to carp_txt using meshtool.\n"
                             "paraview       Opens each file in paraview. Useful for selecting\n"
                             "               stimuli, blocks or scars. (not yet implemented)\n"
                             "mat2ply        Reads .mat file and converts to .ply for manipulation\n"
                             "               in MeshLab (i.e. surface construction)\n"
                        )
    parser.add_argument('--ki',
                        help="Add this optional argument to keep intermediately\n"
                             "created files during tetrahedralisation process."
                             "of the tetrahedralized mesh.", nargs='?',
                        default=False, const=True)
    parser.add_argument('--HULK',
                        help="If command is run on HULK or not", nargs='?',
                        default=False, const=True)
    parser.add_argument("--mat",
                        help="Add this argument if input meshes are .mat format. Assumes all mesh directories\n"
                             "are each in their own directory, bundles in a directory called mat_meshes/")
    args = parser.parse_args()
    run(action=args.action, keep_intmed=args.ki, HULK=args.HULK, mat=args.mat)
    # TODO: implement loop for paraview too (selecting scars and stimuli)
