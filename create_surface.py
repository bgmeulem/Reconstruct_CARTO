"""File containing functions to create a refined two-layered mesh."""

import numpy as np
import pyvista as pv
import pandas as pd
import matplotlib.pyplot as plt
import pymesh as pm
from tqdm import tqdm
from mesh_tools import dist, makePyVista, makePyMesh, cleanMesh

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  # six 'fivethirtyeight' themed colors


def pvToPmCells(pyvistafaces):
    """A pyvista cell looks like [npoints ind1 ind2 ...]
    n_points : amount of points in cell
    Not necessarily a triangle, not necessarily co-planar
    this converts a pyvista face to a pymesh face: [in1 ind2 ind3]"""
    f = []  # list of all faces
    i = 0
    while i < len(pyvistafaces):
        n_points = pyvistafaces[i]
        cell = []
        for co in range(n_points):
            cell.append(pyvistafaces[i + co + 1])
        f.append(cell)
        i += n_points + 1

    return np.array(f)


def getDistances(mesh, type="pymesh"):
    """Gets all edge lengths from a PyMesh mesh (used predominantly in homogenizeMesh())"""
    if type == "pymesh":
        mesh_ = makePyVista(mesh)
    elif type == 'pyvista':
        mesh_ = mesh
    elif type != "pyvista":
        print("Invalid type. Must be \"pymesh\" or \"pyvista\"")
        return 1

    edges = pv.DataSetFilters.extract_all_edges(mesh_)
    pmedges = pvToPmCells(edges.extract_cells(range(edges.n_cells)).cells)  # extract edge ind as cells
    distances = []
    for pair in pmedges:
        co1, co2 = edges.points[pair]
        distances.append(dist(co1, co2))
    return distances


def homogenizeMesh(pvmesh, nsteps=1, boxplot=False, plot_mesh=False, verbose=False, return_dist=False,
                   min_edge=400., max_edge=900.):
    """ Iteratively splits long edges and collapses short ones until they are between
    min_edge and max_edge"""

    # Keep in mind that max allowed distance in split_long_edges should be about twice
    # the tolerance in order to get converging behaviour, since short edges become
    # at least twice as big as tolerance. Setting max_edge < 2 * min_edge wastes some
    # computational time, but can be useful in case the input mesh has self-intersections.

    mesh_ = makePyMesh(pvmesh)
    longest_edge = np.max(getDistances(mesh_))  # used to calculate alpha

    def calcAlpha(n_, nsteps_=nsteps, max_edge_=max_edge, longest_edge_=longest_edge):
        return (1 - n_ / nsteps_) * longest_edge_ * np.exp(-3. * n_ / nsteps_) + max_edge_

    def calcTol(n, nsteps=nsteps, min_edge=min_edge):
        return 500. * (1 - n / nsteps) + min_edge  # min edge length drops linearly

    if boxplot:
        dist_range = [getDistances(mesh_)]  # initialize list of distances during refinement procedure
    for n in tqdm(range(1, nsteps+1), desc='        Resizing edges'):  # don't use \t in tqdm descriptions
        # alpha can never be smaller than 2*tol. Keep some margin in-between
        # It will still run otherwise, but tetgen will probably detect self-intersections

        alpha = calcAlpha(n, nsteps, max_edge, longest_edge)
        tol = calcTol(n, nsteps, min_edge)
        if verbose:
            print("\t\n----------------------------------------------------------")
            print("\tIteration {}".format(n + 1))
            # Splitting long edges
            print("\tSplitting long edges")
        splitmesh_, info = pm.split_long_edges(mesh_, max_edge_length=alpha)

        # Collapsing short edges
        if verbose:
            print("\tCollapsing short edges")
            print("\tTolerance = {}".format(tol))
        colmesh_, info = pm.collapse_short_edges(splitmesh_, tol, preserve_feature=True)

        if verbose:
            dist_col = getDistances(colmesh_)
            print("\t{} edges collapsed\n".format(info['num_edge_collapsed']))
            print("\tMean distance: {}".format(np.mean(dist_col)))
            print("\tSigma: {}".format(np.var(dist_col) ** .5))

        if plot_mesh:
            tri_ = makePyVista(mesh_)  # for plotting
            tri_.plot(show_edges=True, below_color='purple', clim=[-4, 0])
            # tri.slice_orthogonal().plot()
        if boxplot:
            dist_range.append(getDistances(colmesh_))
        mesh_ = colmesh_  # update mesh_ for re-iteration

    if boxplot:
        def plotEdgelengths(dist_range, show=True, save=True):
            """Plots a boxplot of the edgelengths at each iteration step
            dist_range should be a 2D array, each array entry containing all the edge lengths at
            an iteration step"""
            plt.style.use('fivethirtyeight')
            medians = [np.median(d) for d in dist_range]
            means = [np.mean(d) for d in dist_range]
            sigmas = [np.sqrt(v) for v in [np.var(d) for d in dist_range]]

            fig, ax = plt.subplots(figsize=(8, 6))
            ax.set_ylabel("Edge length (µm)", size=22)
            ax.set_xlabel("Iteration step", size=22)
            plt.tight_layout(pad=1.8, h_pad=1.8)
            ax.fill((.5, nsteps+1.5, nsteps+1.5, .5), (max_edge, max_edge, min_edge, min_edge), color=colors[0],
                    alpha=.3, label='Desired edgelength\n({} - {} µm)'.format(min_edge, max_edge))
            # plt.axhline(max_edge, color='black', lw=2.5)
            # plt.axhline(min_edge, color='black', lw=2.5)
            polygonx = [.88] + list(range(2, nsteps + 2)) + list(reversed(range(2, nsteps + 2))) + [.88]
            polygony = [m + s for m, s in zip(means, sigmas)] + \
                       list(reversed([m - s for m, s in zip(means, sigmas)]))
            # ax.fill(polygonx, polygony, color=colors[0], alpha=.2, label="mean $\pm$ $\sigma$")
            ax.boxplot(dist_range, labels=range(len(dist_range)), whis=[0, 100],
                       boxprops=dict(color="gray", linewidth=2.5),
                       medianprops=dict(color="gray", linewidth=2.5),
                       whiskerprops=dict(color="gray", linewidth=2.5),
                       capprops=dict(color="gray", linewidth=2.5), zorder=1)
            ax.plot(range(1, nsteps + 2), [calcAlpha(step) for step in range(nsteps+1)], color=colors[1],
                    )  # xaxis needs to be shifted by one for matplotlib bullshit reasons, hence range(1, nsteps+2)
            ax.plot(range(1, nsteps + 2), [calcTol(step) for step in range(nsteps + 1)], color=colors[1],
                    label=r'$\alpha$ and tolerance')
            ax.set_title("Mesh edge lengths during refinement", size=28)
            ax.legend(loc="upper right", prop={'size': 20})
            ax.tick_params(axis='x', labelsize=18)
            ax.tick_params(axis='y', labelsize=18)
            if show:
                plt.show()
            if save:
                fig.savefig("../Plots/EdgeLengthRefinement_{}-{}.png".format(min_edge, max_edge), dpi=300)

        plotEdgelengths(dist_range)
    print("\tCleaning mesh...")
    mesh_ = cleanMesh(makePyVista(mesh_), tol=min_edge / 2, iter=6)

    if return_dist:
        return mesh_, getDistances(makePyMesh(mesh_))
    return mesh_  # return final version of mesh


def writeVTK(tetmesh, data, names, filename="testmesh.vtk"):
    of = open(filename, "w+")

    #### Header
    of.write("# vtk DataFile Version 3.0\n"
             "vtk output\n"
             "ASCII\n"
             "DATASET UNSTRUCTURED_GRID\n"
             "\n")

    #### Points
    points = tetmesh.points
    of.write("POINTS {} float\n".format(len(points)))
    for p in points:
        for i, co in enumerate(p):
            of.write(str(co))
            if i < 2:
                of.write(" ")
            else:  # after z-component
                of.write("\n")

    #### Cells
    cells_ = pvToPmCells(tetmesh.cells)
    # CELLS {n_cells} {amount of numbers needed to represent one cell times amount of cells}
    of.write("CELLS {} {}\n".format(len(cells_), len(cells_) * 5))
    for cell in cells_:
        of.write("4")
        for node in cell:
            of.write(" {}".format(node))
        of.write("\n")

    #### Cell type
    # 10 is vtk code for tetrahedron
    of.write("CELL_TYPES {}\n".format(len(cells_)))
    for _ in cells_:
        of.write("10\n")

    #### Cell data
    if len(data) > 0:
        for attr, name in zip(data, names):
            of.write("CELL_DATA {}\n"
                     "SCALARS {} float\n"
                     "LOOKUP_TABLE default\n".format(len(cells_), name))
            for _ in attr:
                of.write(str(_) + '\n')  # conductivity tag

    #### Fiber direction
    # of.write("VECTORS fiber float\n")
    # for _ in range(len(cells_)):
    #     of.write("1.000000 0.000000 0.000000\n")

    of.close()


def writeSmesh(mesh, filename="testmesh"):
    of = open(filename + '.smesh', 'w+')
    # attr_names = mesh.array_names ??

    # Write points
    of.write("# Part 1 - the node list\n")
    # <# of points> <dimension (3)> <# of attributes> <boundary markers (0 or 1)>
    # npoints, 3D, 1 attribute (= g_region 1107558400), 1 for boundary marker
    attr_names = ["g_region"]  # hard-coded until correct implementation
    of.write("{} 3 {} 1\n".format(mesh.n_points, len(attr_names)))  # header
    pointbar = tqdm(mesh.points, position=0, desc='        Writing points  ')
    for i, pt in enumerate(pointbar):
        of.write("{} {} {} {} ".format(i, *[round(co, 3) for co in pt]))  # index and co
        of.write("1107558400 ")  # hard-coded g_region as example, corresponds to normal cardiac tissue
        # for attr in attr_names:  # attributes
        # of.write(str(mesh[attr][i]) + " ")
        of.write("1\n")  # boundary marker

    # Write faces
    of.write("# Part 2 - the facet list\n")
    # <# of facets> <boundary markers (0 or 1)>
    of.write("{} 1\n".format(mesh.n_faces))  # header

    faces = pvToPmCells(mesh.faces)
    facesbar = tqdm(faces, position=0, desc='        Writing faces  ')
    for f in facesbar:
        # 1 polygon, boundary marker
        of.write("3")
        for v in f:
            of.write(" {}".format(v))
        of.write(" 1\n")
    pointbar.close()
    facesbar.close()

    # Write hole
    of.write("# Part 3 - the hole list\n")
    # <# of holes>
    of.write("1" + '\n')
    # Coordinate of hole (i.e. center of mesh)
    of.write("0 {} {} {}\n".format(*np.array(mesh.points).mean(axis=0)))


def writeMtr(mesh, length=100., filename="testmesh"):
    """Writes mesh sizing function file for TetGen to use
    Each node should have one value that corresponds to the desired edge length"""
    with open("{}.mtr".format(filename), "w+") as of:
        of.write(str(mesh.n_points) + " 1\n")
        for _ in mesh.points:
            of.write(str(float(length)) + '\n')
    of.close()


def run(of='testmesh', meshdir='', verb=False, boxplot=False, plot_mesh=False, refine_steps=7,
        min_edge=400., max_edge=900., returnmesh=False):
    print("\tReading points and connectivities...")
    triangles = pd.read_csv(meshdir + 'TrianglesSection.csv', sep=',')

    def getPoints(meshdir=meshdir):
        """Reads in endo.txt and epi.txt as generated by AddPoints.py"""
        endo = pd.read_csv(meshdir + 'endo.txt', sep=',')
        epi = pd.read_csv(meshdir + 'epi.txt', sep=',')
        # triangles = triangles[triangles["GroupID"] != -1000000]
        outfile = open("{}.txt".format(of), "w+")
        outfile.write("Index,X,Y,Z\n")

        # Reading points
        points_endo = []
        points_epi = []
        # Multiply all point ordinates by 1000
        # carto .mesh is in mm, openCarp is in µm
        for point in endo[['X', 'Y', 'Z']].to_numpy():
            points_endo.append(list(point))
        for point in epi[['X', 'Y', 'Z']].to_numpy():
            points_epi.append(list(point))
        return [points_endo, points_epi]  # ignore mid, too prone to self-intersection + no reason to have this layer

    def getTriangles(triangles_csv=triangles):
        # Same face for endo and epi
        faces = []
        for row in triangles_csv[['Vertex0', 'Vertex1', 'Vertex2']].to_numpy():
            # face looks like <ind0 ind1 ind2>
            face = []
            for vtx in row:
                face.append(vtx)
            faces.append(face)
        return faces

    def joinLayers(all_points, faces):
        """Uses Pymesh to form meshes and PyVista to add them together"""
        print('\tCreating two-layered mesh...')
        meshes = []
        for points in all_points:
            mesh = pm.form_mesh(points, np.array(faces))
            # for i, name in enumerate(attr_names):
            #     mesh.add_attribute(name)
            #     mesh.set_attribute(name, attrs[i])
            meshes.append(mesh)  # .append method also computes outer hull and forces meshes to be closed (for tetgen)

        fullmeshlist = [makePyVista(mesh) for mesh in meshes]
        fullmesh = fullmeshlist[0]
        # pv.save_meshio('Surface_Backup.stl', fullmesh)  # failsafe, to read in meshlab in case this code fails
        for m in fullmeshlist[1:]:  # sum() operand not supported by PyVista for adding meshes in list, only + operator
            fullmesh += m
        return fullmesh

    all_points = getPoints(meshdir)
    faces = getTriangles(triangles)

    # The only attribute name present in carto .mesh (hardcoded for easy use)
    attr_names = ['GroupID']
    attrs = [np.array(triangles[name]) for name in attr_names]

    fullmesh = joinLayers(all_points, faces)
    fullmesh = cleanMesh(fullmesh, tol=min_edge / 2)
    fineMesh, dist_ = homogenizeMesh(fullmesh, nsteps=refine_steps, boxplot=boxplot,
                                     plot_mesh=plot_mesh, verbose=verb, return_dist=True, min_edge=min_edge,
                                     max_edge=max_edge)

    print("\tWriting .smesh and .mtr file")

    writeSmesh(fineMesh, "{}_{}-{}µm".format(of, int(min_edge), int(max_edge)))
    writeMtr(fineMesh, min_edge, "{}_{}-{}µm".format(of, int(min_edge), int(max_edge)))
    print('\tDone.')
    if returnmesh:
        return fineMesh  # pyvista format
    return 0


if __name__ == '__main__':
    run(meshdir='../carto_meshes/OC4/', boxplot=True, plot_mesh=False, refine_steps=10, returnmesh=False)
    # Reading tetgen-generated .node, .face and .ele files can be done as follows:
    # points, faces, tets = pd.read_csv('testmesh.1.node', usecols=[1, 2, 3],
    #                                   delim_whitespace=True,
    #                                   header=0,
    #                                   names=['x', 'y', 'z']), \
    #                       pd.read_csv('testmesh.1.face', usecols=[1, 2, 3],
    #                                   delim_whitespace=True,
    #                                   header=0,
    #                                   names=['vtx1', 'vtx2', 'vtx3']), \
    #                       pd.read_csv('testmesh.1.ele',
    #                                   usecols=[1, 2, 3, 4],
    #                                   delim_whitespace=True,
    #                                   header=0,
    #                                   names=['Ind', 'vtx1', 'vtx2', 'vtx3', 'vtx4'])
