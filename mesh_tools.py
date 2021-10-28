"""
This file contains a selection of useful vector operations and mesh operations such as io, cleaning, converting etc.
This file is imported by :mod:`carto_mesh`.
"""
import pyvista as pv
import numpy as np
from typing import Tuple, List, Union
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
import time
import matplotlib.pyplot as plt
from collections import OrderedDict
import configparser

plt.style.use('fivethirtyeight')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  # six 'fivethirtyeight' themed colors

G_REGION = 1107558400  # Carto ID for conductive region


def normalize(vector: Union[np.ndarray, List, Tuple]) -> List:
    """
    Normalizes a vector such that its length equals 1.

    Args:
        vector: A 1xN array containing the vector components.

    Returns:
        List: the normalized vector
    """
    return [e / np.sqrt(np.dot(vector, vector)) for e in vector]


def calcAvNormal(facets: pd.DataFrame) -> List:
    """
    Calculates the average normal of a dataframe of triangles. Dataframe must contain
    the components of the normals of each triangle in columns named 'NormalX', 'NormalY', 'NormalZ'
    as is the case for Carto data

    Args:
        facets: A Pandas Dataframe of triangles. Must contain the following columns: 'NormalX', 'NormalY', 'NormalZ',
        'GroupID'

    Returns:
        List: An array containing the components of the normalized average normal of the given triangles.
    """
    av_normal = [0, 0, 0]
    # loop over all facets and calculate average normal
    for row in facets[['NormalX', 'NormalY', 'NormalZ', 'GroupID']].iterrows():
        for i in range(len(av_normal)):  # from 0 to 2
            av_normal[i] += row[1][i]
    av_normal = normalize(av_normal)
    return av_normal


def dist(co1: Union[np.ndarray, List, Tuple], co2: Union[np.ndarray, List, Tuple]) -> float:
    """
    Calculates the distance between two coordinates.

    Args:
        co1: Array containing the components of the first coordinate

        co2: Array containing the components of the second coordinate

    Returns:
        float: Distance between the two coordinates
    """
    return np.linalg.norm([float(e2) - float(e1) for e1, e2 in zip(co1, co2)])


def pmToPvFaces(pmfaces: Union[np.ndarray, List, Tuple]) -> np.ndarray:
    """
    Given an array of pymesh triangles, where these triangles are arrays of the form [ind1, ind2, ind3], converts
    these faces to the PyVista representation of the form [3, ind1_0, ind2_0, ind3_0, 3, ind1_1, ind2_1, ind3_1, 3 ...]

    Args:
        pmfaces: A Nx3 array, where N = n_faces.

    Returns:
        ndarray: A numpy ndarray containing the triangles in PyVista's format
    """
    assert pmfaces.shape[1] == 3, 'At least one cell does not have 3 indices'
    return np.array([[len(f), *f] for f in pmfaces]).flatten()


def pvToPmCells(pyvistafaces: Union[np.ndarray, List, Tuple]) -> np.ndarray:
    """
    Given an array of pyvista faces, converts these faces to the Pymesh representation.

    Args:
        pyvistafaces: An array containing the faces in PyVista format

    Returns:
        ndarray: A numpy ndarray containing the triangles in PyMesh format
    """
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


def pvToPmFaces(pyvistafaces: Union[np.ndarray, List, Tuple]) -> np.ndarray:
    """
    Given an array of pyvista triangles, converts these faces to the PyMesh representation.

    Args:
        pyvistafaces: An array containing the faces in PyVista format.

    Returns:
        ndarray: A numpy ndarray containing the triangles in PyMesh format
    """
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


def makePyVista(pmmesh: pm.Mesh) -> pv.PolyData:
    """
    Converts PyMesh mesh to PyVista mesh (PolyData)

    Args:
        pmmesh: the mesh in PyMesh.Mesh format

    Returns:
        PolyData: The mesh in PyVista's PolyData format
    """
    names = pmmesh.get_attribute_names()
    mesh = pv.PolyData(pmmesh.vertices, pmToPvFaces(pmmesh.faces))
    for a in names:
        mesh[a] = pmmesh.get_attribute(a)
    return mesh


def makePyMesh(pvmesh: pv.PolyData) -> pm.Mesh:
    """
    Converts PyVista mesh to PyMesh mesh

    Args:
        pvmesh: The mesh in PyVista's PolyData format

    Returns:
        Mesh: The mesh as a PyMesh Mesh object
    """
    names = pvmesh.array_names
    mesh = pm.form_mesh(pvmesh.points, pvToPmCells(pvmesh.faces))
    for name in names:
        mesh.add_attribute(name)
        mesh.set_attribute(name, pvmesh[name])
    return mesh


def colorFromCsv(meshdir: str = "") -> pv.PolyData:
    """
    Makes mesh from VerticesSection.csv and TrianglesSection.csv
    Assigns tags in .mesh file as scalar data on this mesh
    Returns mesh with this scalar data.

    Args:
        meshdir: Directory containing VerticesSection.csv and TrianglesSection.csv

    Returns:
        The mesh with scalar data applied in PyVista's PolyData format
    """
    triangles = pd.read_csv(meshdir + 'TrianglesSection.csv', sep=',')
    tri = np.array(triangles[['Vertex0', 'Vertex1', 'Vertex2']].values)

    vertices = pd.read_csv(meshdir + 'VerticesSection.csv', sep=',')
    vert = np.array([[1000. * co for co in p] for p in vertices[['X', 'Y', 'Z']].to_numpy()])

    color = [e for e in triangles['GroupID'].to_numpy()]

    mesh = pv.PolyData(vert, pmToPvFaces(tri))
    mesh['color'] = color
    mesh = mesh.ctp()  # triangle data to point data
    return mesh


def makeNonCollinear(pvmesh: pv.PolyData, edges: Union[np.ndarray, List, Tuple]) -> pv.PolyData:
    """
    Iterates mesh and replades points in the middle of two colinear edges.
    Edges must be 2D array

    Args:
        pvmesh: The input mesh containing colinearities in PyVista's PolyData format

        edges: A list of point index pairs. Each index pair defines a colinear edge.

    Returns:
        PolyData: The mesh with its colinear edges lifted.
    """
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


def cleanMesh_Meshtool(meshname: str, threshold: float = .3, ifmt: str = 'carp_txt', ofmt: str = 'carp_txt') -> None:
    """
    Calls upon meshtool in a shell command to clean a mesh.

    Args:
        meshname: Name of the mesh to be cleaned
        threshold: Quality threshold. Cleaning will continue until this quality is reached, or the max amount of iterations has been reached

        ifmt: The format of the input mesh. May be; carp_txt, carp_bin, vtk, vtk_bin, mmg, neu, purk, obj, off, gmsh,
        stellar, vcflow

        ofmt: The format of the output mesh. May be: carp_txt, carp_bin, vtk, vtk_bin, vtk_polydata, mmg, neu, obj, off,
        stellar, vcflow, ensight_txt

    Returns:
        None: Nothing. Writes out the cleaned mesh in the desired format.
    """
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


def getEdgeLengths(mesh: pm.Mesh) -> List:
    """
    Gets all edge lengths from a PyMesh mesh (used in :func:`~carto_mesh.CartoMesh.homogenizeMesh`)

    Args:
        mesh: The input mesh in PyMesh's Mesh format

    Returns:
        List: A list containing all the lengths of the edges of the input mesh.
    """
    edges = mesh.extract_all_edges()
    pmedges = pvToPmCells(edges.extract_cells(range(edges.n_cells)).cells)  # extract edge ind as cells
    distances = []
    for pair in pmedges:
        co1, co2 = edges.points[pair]
        distances.append(dist(co1, co2))
    return distances


def convertMesh_Meshtool(meshname: str, ifmt: str = 'vtk', ofmt: str = 'carp_txt') -> None:
    """
    Calls upon meshtool in a shell command to convert a mesh.

    Args:
        meshname: Name of the mesh, excluding file extension.

        ifmt: File format of the input mesh. Can be: carp_txt, carp_bin, vtk, vtk_bin, mmg, neu, purk, obj, off, gmsh,
        stellar, vcflow

        ofmt: File format of the output mesh. Can be: carp_txt, carp_bin, vtk, vtk_bin, vtk_polydata, mmg, neu, obj,
        off, stellar, vcflow, ensight_txt

    Returns:
        None: Nothing. Writes out the mesh in the desired format.
    """
    imsh = omsh = meshname
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


def ptsToParaview(filename: str, column_name: str = "meshID") -> None:
    """
    Takes a pts file, adds point ID as extra data to each point, writes out to csv file. This function is useful
    to convert meshes and read them in in ParaView. ParaView does not remember original mesh indexing during mesh
    operations.

    Args:
      filename: Name of .pts file (including '.pts' extension) to be converted to paraview-friendly format

      column_name: Name of the column to contain the mesh indices. Default = "meshID"

    Returns:
        None: Nothing. Writes out the mesh as a .csv file containing the columns 'X', 'Y', 'Z' and <column_name>
    """

    ofname = filename.split(".")[0]
    outfile = open(ofname + "_paraview.csv", "w+")
    df = open(filename)

    print('\n\tWriting {}_paraview.csv'.format(ofname))
    # csv header
    outfile.write("X,Y,Z,{}\n".format(column_name))
    i = 0
    for line in tqdm(df.readlines()[1:], desc='        '):
        for e in line[:-2].split():  # x, y or z co-ordinate
            outfile.write(str(e) + ',')  # write co-ordinates
        outfile.write(str(i) + '\n')  # add point ID
        i += 1

    outfile.close()
    df.close()


def writeToSmesh(mesh: pv.PolyData, name: str) -> None:
    """
    Writes out a .smesh file with a hole in the middle. This is used in :func:`~carto_mesh.CartoMesh.reconstruct` to be
    used as input for TetGen.

    Args:
        mesh: The surface mesh, ready to be tetrahedralized, in PyVista's PolyData format.

        name: Name of the input mesh, and corresponding output .smesh file

    Returns:
        None: Nothing. Writes out a .smesh file.
    """
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


def cleanMesh(pvmesh: pv.PolyData, tol: float, iterations: int = 10, print_si: bool = True) -> pv.PolyData:
    """
    Uses built-in PyVista and PyMesh methods to clean the mesh.

    Args:
        pvmesh: The mesh in PyVista's PolyData format

        tol: Absolute tolerance to use in PyMesh's remove_duplicated_vertices() and PyVista's clean() methods

        iterations: Amount of iterations to try and remove self-intersections with PyMesh

        print_si: Print progress of self-intersection cleanup

    Returns:
        PolyData: The cleaned mesh
    """
    mesh_ = pvmesh.clean(lines_to_points=True, polys_to_lines=True, tolerance=tol, absolute=True)
    mesh_ = makePyMesh(mesh_)
    mesh_, info = pm.remove_degenerated_triangles(mesh_)
    si = len(pm.detect_self_intersection(mesh_))
    if print_si:
        print('\t' + str(si) + ' self-intersections detected')

    i = 0
    while i < iterations and si != 0:  # there are still self-intersections, max 10 iterations
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
    mesh_ = makePyVista(mesh_).clean(lines_to_points=True, polys_to_lines=True, tolerance=tol, absolute=True)
    return mesh_


def getGroupIds(csvfile: str = 'TrianglesSection.csv', col_name: str = "GroupID", skip: Tuple = (0, -1000000)) -> List:
    """
    Extract the tags from TrianglesSection.csv as generated by :func:`~carto_mesh.CartoMesh.__cartoToCsv`

    Args:
        csvfile: Name of the .csv file to read. Default = 'TrianglesSection'. File must be .csv and contain a column <col_name>.

        skip: Array of tags to ignore. These won't be extracted. Default = [0, -1000000]

    Returns:
        List: a list containing the unique tags in the input file
    """
    if skip == "default":
        skip = [0, -1000000]  # skip border region and myocardium
    f = pd.read_csv(csvfile)
    ids = set([e for e in f[col_name].values if e not in skip])
    return list(ids)
