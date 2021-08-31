"""File containing functions to refine tetrahedral mesh using TetGen
Currently unused in reconstruct_mesh.py"""
import pyvista as pv
from subprocess import Popen, PIPE
import argparse
import glob


def writeMtr(meshname, edge_length=None):
    if edge_length is None:
        edge_length = [300.]  # same edge length for all points
    mesh = pv.read(meshname)
    if len(edge_length) == 1:  # same length for all points
        meshname_ = meshname.split('.')
        iteration = int(meshname_[1])  # iteration 1: first output of tetgen, no extra refinement yet.
        with open(meshname_[0] + '.{}.mtr'.format(str(iteration)), 'w+') as of:
            # n_points and metric (always 1)
            of.write(str(mesh.n_points) + " 1\n")
            for _ in range(mesh.n_points):
                of.write(str(float(edge_length[0])) + "\n")
    else:
        assert (len(edge_length) == len(mesh.points)), \
            "Length of edge length array ({}) must be equal to amount " \
            "of mesh points ({}).".format(len(edge_length), mesh.n_points)
        with open('.'.join(meshname.split('.')[:-1]), 'w+') as of:
            for l_ in edge_length:
                of.write(str(int(l_)) + "\n")

    of.close()


def refineMesh(meshname, switches="-rYkANEFq2.5/20.0 -a5e+7"):
    """Refines tetrahedral mesh with TetGen"""
    command = "tetgen {} {}".format(switches, meshname)
    print("\n\tRunning bash command: " + command)
    print("\t-------------------- TetGen output --------------------")
    p = Popen(['stdbuf', '-o0'] + command.split(" "), stdout=PIPE, stderr=PIPE, encoding='utf-8')
    for line in p.stdout:
        print("\t" + str(line.rstrip()))
        p.stdout.flush()
    print("\t------------------ End TetGen output ------------------")


def run(meshname, switches="-rYkANEFq2.5/20.0 -a5e+7", edge_length=None):
    """Re-reads existing .smesh file and refines with TetGen"""
    if edge_length is None:
        edge_length = [300.]
    writeMtr(meshname, edge_length)
    refineMesh(meshname, switches)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--filename',
                        help="Input smesh file, including .smesh format. Default is first .smesh file in cwd.",
                        type=str, default=None)
    parser.add_argument("--switches",
                        help="TetGen switches",
                        type=str, default="-rmYkq1.8/20 -a4e5")
    args = parser.parse_args()
    if not args.filename:
        fn = glob.glob("*.vtk")[0]  # take first .vtk file
    else:
        fn = args.filename
    run(fn, switches=args.switches, edge_length=[300.])
