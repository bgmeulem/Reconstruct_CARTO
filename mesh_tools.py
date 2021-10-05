import random
import pyvista as pv
import numpy as np
from sklearn import neighbors as nb
import pandas as pd
from subprocess import Popen, PIPE
import pymesh as pm
import sys
import io
from tqdm import tqdm


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
    # names = pmmesh.get_attribute_names()
    mesh = pv.PolyData(pmmesh.vertices, pmToPvFaces(pmmesh.faces))
    # for a in names:
    #     mesh[a] = pmmesh.get_attribute(a)
    return mesh


def makePyMesh(pvmesh):
    """Converts PyVista mesh to PyMesh mesh"""
    # names = pvmesh.array_names
    mesh = pm.form_mesh(pvmesh.points, pvToPmCells(pvmesh.faces))
    # for name in names:
    #     mesh.add_attribute(name)
    #     mesh.set_attribute(name, pvmesh[name])
    return mesh


def colorFromCsv(meshdir=""):
    """Makes mesh from VerticesSection.csv and TrianglesSection.csv (generated by carto2csv.py)
    Assigns tags in .mesh file as scalar data on this mesh
    Returns mesh with this scalar data"""
    triangles = pd.read_csv(meshdir+'TrianglesSection.csv', sep=',')
    tri = np.array(triangles[['Vertex0', 'Vertex1', 'Vertex2']].values)

    vertices = pd.read_csv(meshdir+'VerticesSection.csv', sep=',')
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
                    print('\t'+val)
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
                    print('\t'+val)
                    output_map.pop(key)
                    break


def vtkToStl(vtk_mesh, location, meshname):
    # mesh = pv.PolyData(mesh.points, mesh.cells)
    faces = pv.DataSetFilters.extract_surface(vtk_mesh)
    print(faces)
    print(faces.faces)
    pv.save_meshio(location+'/'+meshname.split('.')[0]+'.stl', faces)


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
            "1107558400 ")  # hard-coded g_region as example, corresponds to normal cardiac tissue TODO: magic number
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
            sys.stdout.write("\r" + '\tIteration {}: '.format(i+1) + str(si) + ' self-intersection(s) left ')
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

