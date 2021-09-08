import pyvista as pv
import pandas as pd
import numpy as np
from create_surface import writeVTK, pvToPmCells
import carto2csv
from tqdm import tqdm
import argparse
import matplotlib.pyplot as plt
import glob
import os
from sklearn import neighbors as nb
import colors
from mesh_tools import getGroupIds

plt.style.use("fivethirtyeight")


def writeAdjust(meshdir, mesh):
    """Writes adjustment file to close off Na2+ channels in cells where CV ~ 0"""
    # mesh needs to be pyvista PolyData
    # meshname = glob.glob(meshdir+"*.1.vtk")[0].split('.')[0]
    cells = pvToPmCells(mesh.cells)  # mesh cells in PyMesh format
    speed = mesh["speed"]

    ptn = np.ones(len(mesh.points)) * 7.8  # default value for g_Na
    for i in range(len(speed)):
        if speed[i] < 1e-5:  # if CV is ~ 0 -> close off Na channel
            vertices = cells[i]
            ptn[vertices] = 0.

    stimfile = open(meshdir + "gNA2.adj", "w+")
    stimfile.write(str(len(ptn)) + "\n" + "extra\n")
    for i in range(len(ptn)):
        stimfile.write(str(i) + " " + str(ptn[i]) + "\n")
    stimfile.close()

    datfile = open(meshdir + "gNA2.dat", "w+")
    dat = np.zeros(len(mesh.points))
    dat[ptn == 0.] = 1
    for e in dat:
        datfile.write(str(e) + '\n')
    datfile.close()


def getNonCondReg(meshdir, mesh, plot=False):
    scars = pd.DataFrame(columns=["meshID", "speed", "x_um", "y_um", "z_um"])  # LPV, RPV and MV
    if os.path.isdir(meshdir + "Scars"):  # apply scars if directory exists
        for csvfile in glob.glob(meshdir + "Scars/*.csv"):
            scar = pd.read_csv(csvfile)
            # write scar .dat file
            datfile = open(csvfile.split(".")[0] + ".dat", "w+")
            dat = np.zeros(len(mesh.points))
            dat[scar["meshID"]] = 1
            for e in dat:
                datfile.write(str(e) + '\n')
            datfile.close()
            scars = scars.append(pd.DataFrame([[index, 0., *mesh.points[index]] for index in scar["meshID"].values],
                                              columns=["meshID", "speed", "x_um", "y_um", "z_um"]))
            # plot scars one by one:
            if plot:
                testmesh = mesh
                testmesh["speed"] = mesh.n_points * [0.]
                testmesh["speed"][scar["meshID"]] = 1.
                testmesh.ptc()
                testmesh.plot(stitle=csvfile)
    return scars


def getCartoScar(mesh, scarfile, k=10):
    carto_scar = open(scarfile, "r").readlines()
    id1, scar_co = [line.split(" : ")[0] for line in carto_scar], \
                   [line.split(" : ")[1].strip('\n') for line in carto_scar]
    scar_co = pv.PolyData([[1000. * float(c) for c in p.split(" ")] for p in scar_co])  # mm to µm conversion

    tree = nb.KDTree(mesh.points)
    distances, indices = tree.query(scar_co.points, k)  # k closest mesh points to scar points, i.e. mesh projection
    proj_scar_indices = []  # indices of points on mesh closest to carto scar
    for index in indices.flatten():
        if index not in proj_scar_indices:
            proj_scar_indices.append(index)
    carto_scar = pv.PolyData(mesh.points[proj_scar_indices])
    carto_scar["speed"] = len(carto_scar.points) * [0.]
    return carto_scar, scar_co, proj_scar_indices


def plotMesh(pvdata, mesh, scar_co=None, proj_scar_indices=None):
    plotter = pv.Plotter(title="Conduction velocities")
    cells = pvToPmCells(mesh.cells)
    cell_center = mesh.points[cells].mean(1)
    cx, cy, cz = mesh.points.mean(axis=0)  # center of entire mesh
    mask = cell_center[:, 2] > cz  # everything above the y-plane
    cell_ind = mask.nonzero()[0]
    subgrid = mesh.extract_cells([cell_ind])
    plotter.add_mesh(subgrid, opacity=.5)

    cx, cy, cz = pvdata.points.mean(axis=0)
    mask = pvdata.points[:, 2] > cz
    cell_ind = mask.nonzero()[0]
    subgrid2 = pv.PolyData(pvdata.points[cell_ind])
    subgrid2['speed'] = pvdata['speed'][cell_ind]
    plotter.add_mesh(subgrid2, point_size=6)

    if len(scar_co) > 1 and len(proj_scar_indices) > 1:
        plotter.add_mesh(scar_co, color='red', point_size=7)
        plotter.add_mesh(mesh.points[proj_scar_indices], color='yellow')

    plotter.show()


def run(meshdir="", meshname=None, speed_limit=None, plot_mesh=False, radius=4000, sharpness=1.5, write_csv=False,
        write_VTK_file=False, write_txt=True, plot_dist=False, write_dat=False,
        write_xyz=False, write_adj=False, apply_carto_scar=False, n_variations=10, n_neighbors=5,
        outdir='scale_factors/', plot_scar=False, manual_scar=False):
    # Drop out nonsensical speed limits
    if speed_limit is None:
        speed_limit = [0., 1.4]  # default, in mm/ms

    if meshname is None:
        meshnames = glob.glob(meshdir + "*µm.1.vtk")
        print("\tmeshnames: ", meshnames)
        meshnames.sort()
        meshname = meshnames[0].split("/")[-1]  # take first vtk file
    print("\tMeshname: ", meshname)

    # read in data
    input_data = pd.read_csv(glob.glob(meshdir + "speed.csv")[0],
                             usecols=["speed", "x", "y", "z"])  # "output/movie/output/movie/movie_*.csv"
    med_speed = np.mean(input_data["speed"])
    input_data["x_um"] = [1000. * e for e in input_data["x"]]
    input_data["y_um"] = [1000. * e for e in input_data["y"]]
    input_data["z_um"] = [1000. * e for e in input_data["z"]]
    # create new dataframe for mutable purposes
    calculated_data = pd.DataFrame()
    calculated_data["x_um"] = [1000. * e for e in input_data["x"]]
    calculated_data["y_um"] = [1000. * e for e in input_data["y"]]
    calculated_data["z_um"] = [1000. * e for e in input_data["z"]]
    calculated_data["speed"] = input_data["speed"]
    mesh = pv.read(meshdir + meshname)

    speeds = input_data["speed"]
    for i in range(len(speeds)):
        s = speeds[i]
        if s < speed_limit[0]:
            speeds[i] = speed_limit[0]
        if s > speed_limit[1]:
            speeds[i] = speed_limit[1]
    input_data["speed"] = speeds

    ids = getGroupIds()  # Carto tags (excluding 0 and -10000): correspond to MV, LPV and RPV. At max 3 different tags
    print("\tDetected tags (non-conductive): ", ids)
    for n in range(n_variations):
        if n != 0:  # variations of conduction velocities
            print("\n\t#### Variation ", n)
            # tweak conduction velocities of input CV file
            points = input_data[["x_um", "y_um", "z_um"]].values
            tree = nb.KDTree(points)
            speeds = np.zeros(len(input_data))
            for i in tqdm(range(len(points)), desc='        Calculating new velocities'):
                p = points[i]
                dist, neighbors = tree.query([p], k=n_neighbors)
                neighborCVs = input_data.loc[[int(e) for e in neighbors[0]]]["speed"]
                mean, sigma = np.mean(neighborCVs), np.std(neighborCVs)
                new_cv = np.random.normal(mean, sigma, 1)
                speeds[i] = np.abs(new_cv)
            calculated_data["speed"] = speeds

        colors.run()  # creates colors.csv
        if os.path.exists(meshdir + "noncond.csv"):
            non_cond = np.array(pd.read_csv('noncond.csv').values)
            non_cond = pv.PolyData(non_cond)
            non_cond["speed"] = non_cond.n_points * [0.]
            myo = pv.PolyData(np.array(pd.read_csv('myo.csv').values))
            myo["speed"] = myo.n_points * [1.]
            c = non_cond + myo  # a mask that is 1 where the carto point is myocardium
            # for each meshpoint, find closest carto point
            tree = nb.KDTree(c.points)
            distances, indices = tree.query(mesh.points, k=1)  # these are c indices for each meshpoint
            nc_mesh_ind = [ind for ind in range(len(indices)) if c["speed"][indices[ind]] == 0.]

        else:
            nc_mesh_ind = []
            if not manual_scar:
                print("\n\t!!! No regions detected: manual input needed !!!\n")

        # applying speed limit
        speeds = calculated_data["speed"]
        for i in range(len(speeds)):
            s = speeds[i]
            if s < speed_limit[0]:
                speeds[i] = speed_limit[0]
            if s > speed_limit[1]:
                speeds[i] = speed_limit[1]

        # Create PolyData to use in interpolation
        data = calculated_data
        data["speed"] = speeds

        pvdata = pv.PolyData(np.array([data["x_um"], data["y_um"], data["z_um"]]).T)
        pvdata["speed"] = data["speed"]
        if apply_carto_scar:
            print("\tApplying carto scar")
            carto_scar, scar_co, proj_scar_indices = getCartoScar(mesh,
                                                                  glob.glob(meshdir + "coordinates_scar_p*.txt")[0],
                                                                  k=10)
            pvdata += carto_scar

        # Interpolate on mesh
        print("\tInterpolating on mesh")
        mesh = mesh.interpolate(pvdata, radius=radius, sharpness=sharpness,
                                strategy="null_value", null_value=med_speed,
                                pass_point_arrays=False, pass_cell_arrays=False)

        # Set auto-detected non-conductive regions after interpolation
        mesh["speed"] = [0. if p in nc_mesh_ind else mesh["speed"][p] for p in range(mesh.n_points)]
        # Set manually selected scars to 0 velocity
        if manual_scar:
            scars = getNonCondReg(meshdir, mesh, plot_scar)
            mesh["speed"] = [0. if p in scars["meshID"].values else mesh["speed"][p] for p in range(mesh.n_points)]
        pointdata = mesh["speed"]
        mesh = mesh.ptc()  # point data to cell data
        cell_data = mesh["speed"]

        if plot_mesh:
            mesh.plot()

        # point data: squared speed
        sq_point_data = pd.DataFrame([e ** 2 for e in pointdata], columns=["squared speed"])
        # cell data: squared speed
        sq_cell_data = pd.DataFrame([e ** 2 for e in cell_data], columns=["squared speed"])

        # write to csv file
        if write_csv:
            print("\tWriting squared speed to csv")
            sq_point_data.to_csv(meshdir + outdir + "sq_CV_point_{}.csv".format(n), index_label="PointID")
            sq_cell_data.to_csv(meshdir + outdir + "sq_CV_cell_{}.csv".format(n), index_label="CellID")

        if write_xyz:
            print("\tWriting point cloud")
            of = open(meshdir + outdir + "input_points_{}.txt".format(n), "w+")
            of.write("X,Y,Z,speed\n")
            for i in tqdm(range(len(pvdata.points))):
                p = pvdata.points[i]
                d = pvdata["speed"][i]
                for c in p:
                    of.write(str(c) + ",")
                of.write(str(d) + '\n')
            of.close()

        # write text file to read in during simulation
        if write_txt:
            print("\tWriting txt file: {}scale_factor_{}.txt".format(outdir, n))
            of = open(meshdir + outdir + "scale_factor_{}.txt".format(n), "w+")
            for e in sq_cell_data["squared speed"]:
                of.write(str(e) + '\n')
            of.close()

        # .dat file for visualisation purposes
        if write_dat:
            print("\tWriting dat file: {}scale_factor_{}.dat".format(outdir, n))
            of = open(meshdir + outdir + "scale_factor_{}.dat".format(n), "w+")
            for e in sq_point_data["squared speed"]:
                of.write(str(e) + '\n')
            of.close()

        if write_adj:
            print("\tWriting file: gNA2.adj")
            writeAdjust(meshdir, mesh)

        # write to vtk for inspection in paraview
        if write_VTK_file and n == 0:
            print("\tWriting mesh to {}_CV{}.vtk".format(meshname.split('.')[0], n))
            polyd = pv.UnstructuredGrid(mesh.cells, np.array(len(pvToPmCells(mesh.cells)) * [10]), mesh.points)
            polyd["speed"] = cell_data
            writeVTK(mesh, [mesh["speed"]], ["speed"], meshdir + "{}_CV{}.vtk".format(meshname.split('.')[0], n))

        # plot distribution as violin plot
        if plot_dist:
            to_plot = [input_data[(input_data["speed"] > speed_limit[0]) &
                                  (input_data["speed"] < speed_limit[1])]["speed"].values,
                       [e for e in cell_data if e > 0.]]
            fig, ax = plt.subplots(figsize=(8, 6))
            plt.tight_layout(pad=3.)
            ax.set_ylabel("Conduction velocity (mm/ms)")
            ax.set_title("Interpolating conduction velocities on "
                         "mesh\n Radius={} $\mu m$   Sharpness={}".format(radius, sharpness))
            ax.set_xticks([1, 2])
            ax.set_xticklabels(["Input points", "Mesh cells"])
            plt.violinplot(to_plot, showmedians=True, showextrema=False)
            plt.show()
            fig.savefig("Plots/InterpolationR={}_S={}.png".format(radius, ''.join(str(sharpness).split('.'))), dpi=300)

        if len(ids) < 3 and not manual_scar:
            # don't calculate variations if you have to recalculate again with manual scars
            print("\n\t!!! Only {} out of 3 tag(s) found - Manual input needed !!!\n".format(len(ids)))
            break
    to_delete = ["TrianglesSection.csv", "VerticesSection.csv", 'VerticesAttributesSection.csv',
                 'VerticesColorsSection.csv']
    for trash in to_delete:
        if glob.glob(trash):
            os.remove(trash)
            print("\tDeleted ", trash)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--meshdir',
                        help="Mesh directory where .vtk file is, starting from cwd.",
                        type=str, default="")
    parser.add_argument('--outdir',
                        help='Directory where scale factor files will be written to.',
                        type=str, default='scale_factors/')
    parser.add_argument('--msh',
                        help="Name of .vtk file (including .vtk extension)",
                        type=str, default="")
    parser.add_argument('--p',
                        help="plot mesh with interpolated velocity values and sampling nodes.",
                        nargs="?", default=False, const=True)
    parser.add_argument('--nvar',
                        help="Amount of different conduction velocity distributions. First one is always\n"
                             "the unchanged input velocity distribution.",
                        type=int, default=10)
    parser.add_argument('--nnb',
                        help="Amount of neighboring points to infer conduction velocity from.",
                        type=int, default=10)
    parser.add_argument('--radius',
                        help='interpolation radius',
                        type=float, default=5000.)
    parser.add_argument('--sharpness',
                        help='Sharpness of gaussian interpolation distribution',
                        type=float, default=.8)
    parser.add_argument('--writeVTK',
                        help='Write out VTK with squared conductivity as cell data.',
                        nargs='?', default=False, const=True)
    parser.add_argument('--plotscar',
                        help='Write out VTK with squared conductivity as cell data.',
                        nargs='?', default=False, const=True)
    parser.add_argument('--ms',
                        help='Apply manually selected non-conductive regions.\n'
                             'Input = .csv files in \'Scars/\' directory)',
                        nargs='?', default=False, const=True)
    parser.add_argument('--cartoscar',
                        help="apply carto scar. input = \'coordinates_scar_p*.txt\'",
                        nargs='?', default=False, const=True)
    args = parser.parse_args()
    if not args.msh:
        mesh = glob.glob("*µm.1.vtk")[0]
    else:
        mesh = args.mesh

    # create csv files:
    print("\n\tMaking .csv files\n")
    mn = glob.glob("*.mesh")[0]
    carto2csv.cartoToCsv(mn, "")

    # TODO: BEWARE: writeVTK is set to True by hardcode default for current project
    run(meshdir=args.meshdir, meshname=mesh, write_txt=True, write_dat=True, n_variations=args.nvar,
        n_neighbors=args.nnb, write_VTK_file=True, apply_carto_scar=args.cartoscar, radius=args.radius,
        sharpness=args.sharpness, outdir=args.outdir, plot_scar=args.plotscar, manual_scar=args.ms)
