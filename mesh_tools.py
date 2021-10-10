import pyvista as pv
import numpy as np
from sklearn import neighbors as nb
import pandas as pd
from subprocess import Popen, PIPE
import pymesh as pm
import csv
import glob
import os
import sys
import io
from tqdm import tqdm

G_REGION = 1107558400  # Carto ID for conductive region


def normalize(vector):
    return [e / np.sqrt(np.dot(vector, vector)) for e in vector]


def calcAvNormal(facets):
    """Calculates the average normal of a dataframe of triangles. Dataframe must contain
    the components of the normals of each triangle in columns named 'NormalX', 'NormalY', 'NormalZ'
    as is the case for Carto data"""
    av_normal = [0, 0, 0]
    # loop over all facets and calculate average normal
    for row in facets[['NormalX', 'NormalY', 'NormalZ', 'GroupID']].iterrows():
        for i in range(len(av_normal)):  # from 0 to 2
            av_normal[i] += row[1][i]
    av_normal = normalize(av_normal)
    return av_normal


def meshToCsv(inmesh, meshdir="", verbose=True):
    with io.open(inmesh, 'r', encoding='utf-8') as f:
        # variables used throughout loop:
        section = "Header"  # name of current section
        current_line = 0  # index of current line being read
        section_line = 0  # index of current section name
        section_ind = 0

        # Section 0 = Header, file format info
        # Section 1 = General Attributes, info about section data
        # Section 2 = VerticesSection, x y z coordinates plus data (normals and GroupID)
        # Section 3 = TrianglesSection,
        # section 4 = VerticesColorSection
        # Section 5 = VerticesAttributesSection

        for line in f:
            current_line += 1  # update line index

            if line[0] == "[":  # new section encountered
                if verbose:
                    print("\t" + str(section_ind) + '. ' + section)
                section = line[1:-2]  # section name
                section_line = current_line
                section_ind += 1

                # Opening or closing csv files
                if section_ind > 2:
                    of.close()  # past VerticesSection, an outfile has already been created (infra)
                if section_ind > 1:  # skip GeneralAttributes
                    of = open(meshdir + section + ".csv", "w+", encoding='utf-8')

                # New tqdm loops per section
                # if section_ind > 1:
                #     if section_ind > 2:
                #         time.sleep(1e-6)
                #         loop.close()  # close previous loop
                #     if section_ind == 3:
                #         time.sleep(1e-6)
                #         loop = tqdm(desc='        ' + section, total=n_triangles, position=0, leave=True)
                #     else:
                #         time.sleep(1e-6)
                #         loop = tqdm(desc='        ' + section, total=n_vertices, position=0, leave=True)

            # extract n_vertices and n_triangles to use for tqdm loops
            # elif section_ind == 1:
            #     if current_line == section_line + 3:  # n_verticess
            #         n_vertices = int(line.split("=")[1])
            #     elif current_line == section_line + 4:
            #         n_triangles = int(line.split('=')[1])

            elif section_ind > 1:  # useful data
                if section_ind == 4:
                    header_line = section_line + 2  # VerticesColorSection has extra line
                else:
                    header_line = section_line + 1

                if current_line == header_line:  # column names
                    column_names = line.split()[1:]  # first element is useless ";"
                    of.write("Index")
                    for name in column_names:
                        of.write("," + name)
                    of.write("\n")
                    # time.sleep(1e-8)

                elif len(line) > 1 and line[0] != ";":  # actual data, not empty line or column names
                    # time.sleep(1e-8)
                    # loop.update(1)
                    ind, data = line.split("=")  # line has shape "ind = x  y  z  nx  ny  nz  groupID
                    of.write(str(ind))
                    for el in data.split():  # ignore ind
                        of.write("," + str(el))
                    of.write("\n")
        # time.sleep(1e-8)
        # loop.close()
        # time.sleep(1e-8)
        of.close()


def dist(co1, co2):
    return np.linalg.norm([float(e2) - float(e1) for e1, e2 in zip(co1, co2)])


def pmToPvFaces(pmfaces):
    assert pmfaces.shape[1] == 3, 'At least one cell does not have 3 indices'
    return np.array([[len(f), *f] for f in pmfaces]).flatten()


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


def pvToPmFaces(pyvistafaces):
    """A pyvista cell array looks like [npoints_1 ind1_1 ind2_1 ... npoints_i ind1_i ind2_i ...]
    n_points : amount of points in face
    ONLY works for triangles
    this converts a pyvista face to a pymesh face: [in1 ind2 ind3]
    Identical to pvToPmCells, but with extra check to see if all the faces have 3 nodes"""
    f = []  # list of all faces
    i = 0
    while i < len(pyvistafaces):
        n_points = pyvistafaces[i]
        if n_points == 3:  # ignore all non-triangle surfaces
            face_ = []
            for co in range(n_points):
                face_.append(pyvistafaces[i + co + 1])
            f.append(face_)
            i += n_points + 1
        else:
            i += n_points + 1

    return np.array(f)


def makePyVista(pmmesh):
    """Converts PyMesh mesh to PyVista mesh"""
    names = pmmesh.get_attribute_names()
    mesh = pv.PolyData(pmmesh.vertices, pmToPvFaces(pmmesh.faces))
    for a in names:
        mesh[a] = pmmesh.get_attribute(a)
    return mesh


def makePyMesh(pvmesh):
    """Converts PyVista mesh to PyMesh mesh"""
    names = pvmesh.array_names
    mesh = pm.form_mesh(pvmesh.points, pvToPmCells(pvmesh.faces))
    for name in names:
        mesh.add_attribute(name)
        mesh.set_attribute(name, pvmesh[name])
    return mesh


def colorFromCsv(meshdir=""):
    """Makes mesh from VerticesSection.csv and TrianglesSection.csv (generated by carto2csv.py)
    Assigns tags in .mesh file as scalar data on this mesh
    Returns mesh with this scalar data"""
    triangles = pd.read_csv(meshdir + 'TrianglesSection.csv', sep=',')
    tri = np.array(triangles[['Vertex0', 'Vertex1', 'Vertex2']].values)

    vertices = pd.read_csv(meshdir + 'VerticesSection.csv', sep=',')
    vert = np.array([[1000. * co for co in p] for p in vertices[['X', 'Y', 'Z']].to_numpy()])

    color = [e for e in triangles['GroupID'].to_numpy()]

    mesh = pv.PolyData(vert, pmToPvFaces(tri))
    mesh['color'] = color
    mesh = mesh.ctp()  # triangle data to point data
    return mesh


def makeNonCollinear(pvmesh, edges):
    """Iterates mesh and moves over points in the middle of two colinear edges.
    edges must be 2D array"""
    mesh = pvmesh
    points = mesh.points
    triangles = pd.DataFrame(pvToPmFaces(pvmesh.faces), columns=[['Vertex0', 'Vertex1', 'Vertex2']])
    tree = nb.KDTree(mesh.points)
    for edgepair in edges:
        edge1, edge2 = edgepair[0], edgepair[1]  # indices of edges
        cp = np.intersect1d(edge1, edge2)[0]  # index of common point
        p1, p2 = [e for e in edge1 + edge2 if e != cp]
        co_cp = points[cp]  # coordinate of common point

        ### check which one is in the middle:
        segm = [cp, p1, p2]  # indices of segment complex
        # start out with assumption that common point is in the middle
        # and that the sum of the distances to its neighbors is smallest
        minsum = dist(points[cp], points[p1]) + dist(points[cp], points[p2])
        for p in segm:  # now check if that's true
            other_points = [e for e in segm if e != p]
            sum = np.sum([dist(points[p], points[e]) for e in other_points])
            if sum < minsum:  # this point is in the middle!
                cp = p  # cp now stand for "center point" instead of common point
                # this is the point we need to move

        adjacent_faces = triangles.loc[(triangles[['Vertex0']].values == cp) |
                                       (triangles[['Vertex1']].values == cp) |
                                       (triangles[['Vertex2']].values == cp)]

        av_dir = [0, 0, 0]
        for point in set([e for e in adjacent_faces.to_numpy().flatten() if e != cp]):
            co = points[point]
            direction = [p - c for p, c in zip(co, co_cp)]
            av_dir = [a + n for a, n in zip(av_dir, direction)]
        av_dir = [e / np.linalg.norm(av_dir) for e in av_dir]

        # get closest neighbour to center point of colinear edges
        distance, ind = tree.query([co_cp], k=2)  # check edge-range for an idea of this value
        # Move point just a tad aside
        points[cp] = np.array([p + d * distance[0][1] / 10 for p, d in zip(co_cp, av_dir)])

        # move other points in the other direction (optional):
        for p in [e for e in segm if e != cp]:
            points[p] = np.array([e - d * distance[0][1] / 10 for e, d in zip(points[p], av_dir)])
    newmesh = pv.PolyData(points, pvmesh.faces)
    return newmesh


def cleanMesh_Meshtool(meshname, threshold=.3, ifmt='carp_txt', ofmt='carp_txt'):
    p = Popen("meshtool clean quality -msh={0} -thr={1}"
              "-ifmt={2} -ofmt={3} -outmsh={0}".format(meshname, threshold, ifmt, ofmt),
              shell=True, stdout=PIPE, encoding='utf-8')

    # Output category check
    output_map = {'Reading elements': 'Reading elements in text CARP format...',
                  'Reading points': 'Reading points in text CARP format...',
                  'Reading fibers': 'Reading fibers (1 direction) in text CARP format...',
                  'Writing elements': 'Writing elements in text CARP format...',
                  'Writing points': 'Writing points in text CARP format...',
                  'Writing fibers': 'Writing fibers in text CARP format...',
                  'Reading mesh': 'Reading mesh...',
                  'Writing mesh': 'Writing mesh...',
                  'Iteration': 'Shifting points...'}

    # output handling
    for line in p.stdout:
        line = str(line.rstrip())
        if line:
            for key, val in output_map.items():
                if key in line:
                    print('\t' + val)
                    output_map.pop(key)
                    break


def getEdgeLengths(mesh):
    """Gets all edge lengths from a PyMesh mesh (used in homogenizeMesh())"""
    edges = mesh.extract_all_edges()
    pmedges = pvToPmCells(edges.extract_cells(range(edges.n_cells)).cells)  # extract edge ind as cells
    distances = []
    for pair in pmedges:
        co1, co2 = edges.points[pair]
        distances.append(dist(co1, co2))
    return distances


def convertMesh_Meshtool(meshname, ifmt='vtk', ofmt='carp_txt'):
    imsh = meshname
    omsh = meshname
    if ifmt == 'vtk':
        imsh += '.1'
    if ofmt == 'vtk':
        omsh += '.1'
    p = Popen("meshtool convert -imsh={} -ifmt={} -omsh={} -ofmt={}".format(imsh, ifmt, omsh, ofmt),
              shell=True, stdout=PIPE, stderr=PIPE, encoding='utf-8')

    # checks if meshtool is at a certain point in output
    # Output category check
    output_map = {'Reading elements': 'Reading elements in text CARP format...',
                  'Reading points': 'Reading points in text CARP format...',
                  'Reading fibers': 'Reading fibers (1 direction) in text CARP format...',
                  'Writing elements': 'Writing elements in text CARP format...',
                  'Writing points': 'Writing points in text CARP format...',
                  'Writing fibers': 'Writing fibers in text CARP format...',
                  'Reading vtk': 'Reading vtk file...',
                  'Writing vtk file in text format': 'Writing vtk file in text format...',
                  'Reading mesh': 'Reading mesh...',
                  'Writing mesh': 'Writing mesh...'}

    # output handling
    for line in p.stdout:
        line = str(line.rstrip())
        if line:
            for key, val in output_map.items():
                if key in line:
                    print('\t' + val)
                    output_map.pop(key)
                    break


def ptsToParaview(filename: str):  # .pts meshfile to paraview csv
    """
    Takes a pts file, adds point ID as extra data to each point, writes out to csv file
    Args:
      filename: Name of .pts file (including '.pts' extension) to be converted to paraview-friendly format
    """

    ofname = filename.split(".")[0]
    outfile = open(ofname + "_paraview.csv", "w+")
    df = open(filename)

    print('\n\tWriting {}_paraview.csv'.format(ofname))
    # csv header
    outfile.write("X,Y,Z,meshID\n")
    i = 0
    for line in tqdm(df.readlines()[1:], desc='        '):
        for e in line[:-2].split():  # x, y or z co-ordinate
            outfile.write(str(e) + ',')  # write co-ordinates
        outfile.write(str(i)+'\n')  # add point ID
        i += 1

    outfile.close()
    df.close()
    return 0


def vtkToStl(vtk_mesh, location, meshname):
    # mesh = pv.PolyData(mesh.points, mesh.cells)
    faces = pv.DataSetFilters.extract_surface(vtk_mesh)
    pv.save_meshio(location + '/' + meshname.split('.')[0] + '.stl', faces)


def writeToSmesh(mesh, name):
    of = open(name + '.smesh', 'w+')
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
        of.write(
            f"{G_REGION} ")
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


def cleanMesh(pvmesh, tol, iter=10, print_si=True):
    mesh_ = pvmesh.clean(lines_to_points=True, polys_to_lines=True,
                         tolerance=tol, absolute=True)
    mesh_ = makePyMesh(mesh_)
    mesh_, info = pm.remove_degenerated_triangles(mesh_)
    si = len(pm.detect_self_intersection(mesh_))
    if print_si:
        print('\t' + str(si) + ' self-intersections detected')

    i = 0
    while i < iter and si != 0:  # there are still self-intersections, max 10 iterations
        mesh_, info = pm.remove_duplicated_faces(mesh_)
        mesh_, info = pm.remove_duplicated_vertices(mesh_, tol=tol)
        mesh_ = pm.resolve_self_intersection(mesh_)
        si = len(pm.detect_self_intersection(mesh_))
        if print_si:
            sys.stdout.write("\r" + '\tIteration {}: '.format(i + 1) + str(si) + ' self-intersection(s) left ')
            sys.stdout.flush()
        i += 1
    print("")
    mesh_, info = pm.remove_duplicated_vertices(mesh_, tol=tol)
    mesh_, info = pm.remove_duplicated_faces(mesh_)
    mesh_ = makePyVista(mesh_).clean(lines_to_points=True, polys_to_lines=True,
                                     tolerance=tol, absolute=True)
    return mesh_


def getGroupIds(csvfile='TrianglesSection.csv', skip="default"):
    if skip == "default":
        skip = [0, -1000000]  # skip border region and myocardium
    f = pd.read_csv(csvfile)
    ids = set([e for e in f["GroupID"].values if e not in skip])
    return list(ids)


def applyCV(meshdir="", meshname=None, speed_limit=None, plot_mesh=False, radius=4000, sharpness=1.5, write_csv=False,
            write_VTK_file=False, write_txt=True, write_dat=False,
            write_xyz=False, write_adj=False, apply_carto_scar=False, n_variations=10, n_neighbors=5,
            outdir='scale_factors/', plot_scar=False, manual_scar=False):
    def extractTags(meshdir="", writeVtk=False, write_scar=True):
        """Writes out .csv file with scar coördinates. Can also write out this surface as .vtk file"""

        # get mesh from GroupID's
        def writePerTag(mesh, meshdir):
            """Writes out all points in a mesh with a certain tag to a .csv file
            Currently only makes distinction between myocardium and not myocardium aka tag 0 or not 0"""
            print('Tags: ', sorted(set(mesh['color']))[::-1])
            scar = mesh.points[[mesh['color'] != 0]]
            myo = mesh.points[[mesh['color'] == 0]]
            if len(scar) > 0:
                with open(meshdir + 'noncond.csv', 'w+') as of:
                    csvWriter = csv.writer(of, delimiter=',')
                    csvWriter.writerows(scar)
            with open(meshdir + "myo.csv", 'w+') as of:
                csvWriter = csv.writer(of, delimiter=',')
                csvWriter.writerows(myo)

        mesh = colorFromCsv(meshdir)
        if write_scar:
            # write tagged LPV, RPV, MV to noncond.csv and myocardium to myo.csv
            writePerTag(mesh, meshdir)
        if writeVtk:
            pv.save_meshio('colors.vtk', mesh)
        print("Tagged .vtk file written")

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

        extractTags()  # creates colors.csv
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
        if write_VTK_file:
            print("\tWriting mesh to {}_CV{}.vtk".format(meshname.split('.')[0], n))
            polyd = pv.UnstructuredGrid(mesh.cells, np.array(len(pvToPmCells(mesh.cells)) * [10]), mesh.points)
            polyd["speed"] = cell_data
            pv.save_meshio(meshdir + "{}_CV{}.vtk".format(meshname.split('.')[0], n), mesh)

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
