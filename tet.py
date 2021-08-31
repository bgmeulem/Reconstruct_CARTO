from subprocess import Popen, PIPE
import glob
import argparse


def run(filename="1-AT1230CL", switches="-pYkANEFq2.5/20.0 -a5e+7", HULK=False):
    """Runs TetGen in terminal and captures output."""
    tetdir = "~/.config/tetgen1.6.0/build/"  # directory where tetgen is on HULK
    command = "tetgen {} {}.smesh".format(switches, filename)
    if HULK:
        command = tetdir + command
    print('\tNow running bash command {}\n'.format("\"" + command + "\""))
    # print("\t-p  For tetrahedralization\n"
    #       "\t-Y  for preserving input faces as boundaries\n"
    #       "\t-k  Outputs mesh to .vtk file for viewing by Paraview\n"
    #       "\t-A  Assigns region attributes (i.e. conductivity regions)\n"
    #       "\t-N  Suppresses creation of .node file\n"
    #       "\t-E  Suppresses creation of .ele file\n"
    #       "\t-F  Suppresses creation of .face file\n"
    #       "\t-q  minVolume/minDihAngle - Refines mesh (improves qualiy)\n"
    #       "\t-a  volume - Volume constraint\n"
    #       "\t-V  (optional) for verbose output (Mesh quality statistics, edge lengths, volumes...)")
    # print("\n\tA good mesh should have between 5e5 and 1e6 tetrahedra\n")
    print("\t-------------------- TetGen output --------------------")
    # adding stdbuf -o0 makes TetGen output real-time
    if HULK:
        p = Popen('stdbuf  -o0 ' + command, stdout=PIPE, stderr=PIPE, encoding='utf-8', shell=True)
        colin_pair = []
        for line in p.stdout:
            l_ = str(line.rstrip())
            print("\t" + l_)
            # The following checks assume that, if an error is found, it's colinearity
            # This still needs to be expanded for overlapping lines and segments intersecting triangles
            if '1st' in l_ and 'facet' not in l_ and 'segment' not in l_ and 'seg' not in l_:
                # line is of form "1st: [24578,24581]."
                # split at ":", drop first two (whitespace and "[") and last two ("]" and ".") chars
                # split at comma to make list again
                edge1 = [int(e) for e in l_.split(':')[1][2:-2].split(',')]
            if '2nd' in l_ and 'facet' not in l_ and 'segment' not in l_ and 'seg' not in l_:
                # line is of form "2nd: [24578,24581]."
                edge2 = [int(e) for e in l_.split(':')[1][2:-2].split(',')]
                colin_pair.append([edge1, edge2])
            p.stdout.flush()
    else:
        p = Popen(['stdbuf', '-o0'] + command.split(" "), stdout=PIPE, stderr=PIPE, encoding='utf-8')
        colin_pair = []
        for line in p.stdout:
            l_ = str(line.rstrip())
            print("\t" + l_)
            # The following checks assume that, if an error is found, it's colinearity
            # This still needs to be expanded for overlapping lines and segments intersecting triangles
            if '1st' in l_ and 'facet' not in l_ and 'segment' not in l_ and 'seg' not in l_:
                # line is of form "1st: [24578,24581]."
                # split at ":", drop first two (whitespace and "[") and last two ("]" and ".") chars
                # split at comma to make list again
                edge1 = [int(e) for e in l_.split(':')[1][2:-2].split(',')]
            if '2nd' in l_ and 'facet' not in l_ and 'segment' not in l_ and 'seg' not in l_:
                # line is of form "2nd: [24578,24581]."
                edge2 = [int(e) for e in l_.split(':')[1][2:-2].split(',')]
                colin_pair.append([edge1, edge2])
            p.stdout.flush()
    print("\t------------------ End TetGen output ------------------")
    return colin_pair


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    # fetching the arguments
    parser.add_argument('--filename',
                        help="Input smesh file, including .smesh format. Default is first .smesh file in cwd.",
                        type=str, default=None)
    parser.add_argument("--switches",
                        help="TetGen switches",
                        type=str, default="-pYkVq1.5/20 -a5e6")
    parser.add_argument("--HULK",
                        help="Run on HULK?", nargs='?',
                        default=False, const=True)
    args = parser.parse_args()
    if not args.filename:
        fn = glob.glob("*.smesh")[0].split(".")[0]
    else:
        fn = args.filename
    run(fn, switches=args.switches, HULK=args.HULK)  # run with default commands
