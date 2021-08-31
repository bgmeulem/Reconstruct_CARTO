import random
import pyvista as pv
import numpy as np
from sklearn import neighbors as nb
import pandas as pd
from subprocess import Popen, PIPE
import pymesh as pm
import matplotlib.pyplot as plt
import sys
import io
import glob
import os
from scipy.spatial.transform import Rotation
from collections import OrderedDict
# sys.path.append('/'.join(os.getcwd().split('/')[:-1]))  # so that Mesh_Reconstruction directory can be imported
# from Mesh_Reconstruction.carto2csv import meshToCsv #Mesh_Reconstruction. # import can be difficult when in terminal/IDE -> add/delete parent

plt.style.use('fivethirtyeight')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  # six fivethirtyeight themed colors


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
    return np.array([[len(f), *f] for f in pmfaces])


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
    el = False
    wel = False
    po = False
    wp = False
    f = False
    wf = False

    # output handling
    for line in p.stdout:
        line = str(line.rstrip())
        if line:
            if 'Reading elements' in line:
                if not el:
                    el = True
                    print('\tReading elements in text CARP format...')
            elif 'Reading points' in line:
                if not po:
                    po = True
                    print('\tReading points in text CARP format...')
            elif 'Reading fibers' in line:
                if not f:
                    f = True
                    print('\tReading fibers (1 directions) in text CARP format...')
            elif 'ETA' in line:
                pass
            elif 'Writing elements' in line:
                if not wel:
                    wel = True
                    print('\tWriting elements in text CARP format..')
            elif 'Writing points' in line:
                if not wp:
                    wp = True
                    print('\tWriting points in text CARP format...')
            elif 'Writing fibers' in line:
                if not wf:
                    wf = True
                    print('\tWriting fibers (1 direction) in text CARP format...')
            else:
                print('\t' + line)
            p.stdout.flush()


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
    # for vtk to carp conversion
    vtkReading = False
    el = False
    po = False
    f = False
    # for carp to vtk conversion
    rel = False
    rpo = False
    rfi = False
    wvtk = False
    # output handling
    for line in p.stdout:
        line = str(line.rstrip())
        if line:
            # If statements is to prevent messy output (meshtool convert has lots of progress bars)
            # They reduce duplicate output lines to one line
            if 'Reading vtk' in line:
                if not vtkReading:
                    vtkReading = True
                    print('\tReading vtk file...')
            elif 'Writing elements' in line:
                if not el:
                    el = True
                    print('\tWriting elements in txt CARP format...')
            elif 'Writing points' in line:
                if not po:
                    po = True
                    print('\tWriting points in text CARP format...')
            elif 'Writing fibers' in line:
                if not f:
                    f = True
                    print('\tWriting fibers (1 direction) in text CARP format...')
            elif 'Reading elements in text CARP format' in line:
                if not rel:
                    rel = True
                    print('\tReading elements in text CARP format...')
            elif 'Reading points in text CARP format' in line:
                if not rpo:
                    rpo = True
                    print('\tReading points in text CARP format...')
            elif 'Reading fibers (1 directions) in text CARP format' in line:
                if not rfi:
                    rfi = True
                    print("\tReading fibers (1 directions) in text CARP format")
            elif 'Writing vtk file in text format' in line:
                if not wvtk:
                    wvtk = True
                    print('\tWriting vtk file in text format...')
            elif not 'ETA' in line:
                print('\t' + line)
            p.stdout.flush()


def vtkToStl(vtk_mesh, location, meshname):
    # mesh = pv.PolyData(mesh.points, mesh.cells)
    faces = pv.DataSetFilters.extract_surface(vtk_mesh)
    print(faces)
    print(faces.faces)
    pv.save_meshio(location+'/'+meshname.split('.')[0]+'.stl', faces)


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


def getRegionCenters(case='OC59', plot=False, method='mean', location='../carto_meshes'):
    if case == 'all':  # for testing purposes
        dirs = glob.glob('../carto_meshes/*/')
        meshnames = glob.glob("../carto_meshes/*/*µm.1.vtk")
        meshnames.sort(key=lambda x: int(x.split('/')[2][2:]))
        dirs.sort(key=lambda x: int(x.split('/')[2][2:]))
    else:
        meshnames = glob.glob('{}/{}/*µm.1.vtk'.format(location, case))
        dirs = glob.glob('{}/{}/'.format(location, case))

    for m, d in zip(meshnames, dirs):
        regions = ['MV', 'LPV', 'RPV']
        auto_detected = []
        center_points = {}  # to return

        for name in regions:  # manually selected tags in Scars/ directory
            if os.path.exists(d+'Scars/{}.csv'.format(name)):
                scars = [e for e in glob.glob(d+'Scars/*.csv') if name in e]
                scar = pd.DataFrame(columns=['meshID'])
                for s in scars:
                    scar = scar.append(pd.read_csv(s))
                mesh = pv.read(m)
                if method == 'mean':
                    cp = np.mean(mesh.points[[e[0] for e in scar.values]], axis=0)
                elif method == 'median':
                    cp = np.median(mesh.points[[e[0] for e in scar.values]], axis=0)
                center_points[name] = cp.tolist()
            else:
                auto_detected.append(name)

        if len(auto_detected):  # not all regions were manually selected: some were present and functional in .mesh file
            if case in ['OC52']:
                tagdict = {'MV': [-8], 'LPV': [-5, -4], 'RPV': [-6, -7]}
            else:
                tagdict = {'MV': [-4], 'LPV': [-5], 'RPV': [-6]}
            meshToCsv(glob.glob(d+'*.mesh')[0], meshdir=d, verbose=False)  # make tagged vtk file

            cm = colorFromCsv(d)  # tagged (color) mesh
            for name in auto_detected:  # tags present in CARTO .mesh file
                t = tagdict[name]
                points = cm.points[np.in1d(cm['color'], t)]  # get all points that have tag t
                if method == 'mean':
                    cp = np.mean(points, axis=0)  # center point of one tag
                elif method == 'median':
                    cp = np.median(points, axis=0)
                center_points[name] = cp.tolist()

        if plot:
            p = pv.Plotter()
            p.add_mesh(cm, clim=[-8, 0], below_color='grey')
            p.add_mesh(pv.PolyData(np.array([np.array(e) for e in center_points.values()])), color='red', point_size=20)
            p.show()
        # return as dictionary, sorted by key so it's always in the same order
        return OrderedDict(sorted(center_points.items()))


def getBandsAroundRegions(case='OC59', k=1000, meshname='', location='../carto_meshes'):
    """Returns a dict with as keys the regions and as values the indices of the closest
    myocardial points to these regions
    Indices correspond to input mesh(name)
    Make sure the input mesh is aligned with the .vtk mesh"""
    d = location + '/' + case + '/'
    if meshname:
        if '.csv' in meshname:
            mesh = pv.PolyData(pd.read_csv(meshname)[['X', 'Y', 'Z']].values)
            vtkmesh = pv.read(glob.glob(d + '*µm.1.vtk')[0])  # self-selected scar indices are encoded in vtk mesh
        elif '.vtk' in meshname:
            mesh = pv.read(meshname)
        else:
            assert os.path.isdir(meshname), "Mesh {} not found".format(meshname)
            print("Mesh {} not of .csv or .vtk format")

    regions = ['MV', 'LPV', 'RPV']
    auto_detected = []
    close_points = {}  # to return

    for name in regions:  # manually selected tags in Scars/ directory
        if not os.path.exists(d + 'Scars/{}.csv'.format(name)):
            # If self-selected scar doesn't exist, assume it's because it's auto-detected
            auto_detected.append(name)
        else:  # manually selected scar region
            scars = [e for e in glob.glob(d+'Scars/*.csv') if name in e]
            for s in scars:
                scar = pd.read_csv(s)
                scar_ind = [e for e in range(len(mesh.points)) if e in scar['meshID'].values]
                scar_ind = random.sample(scar_ind, len(scar_ind) // 10)
                if '.vtk' in meshname:
                    ind = [e for e in range(len(mesh.points)) if e not in scar['meshID'].values]
                    cut_mesh = pv.PolyData(mesh.points[ind])  # cut away scar
                    scar_mesh = pv.PolyData(mesh.points[scar_ind])
                    cut_mesh['Id'] = ind
                    tree = nb.KDTree(cut_mesh.points)
                    dist, nbs = tree.query(scar_mesh.points, k=k)
                    nbs = list(set(np.array(nbs).flatten()))
                    nbs = cut_mesh['Id'][nbs]
                elif '.csv' in meshname:
                    scar_mesh = pv.PolyData(vtkmesh.points[scar_ind])
                    # p = pv.Plotter()
                    # p.add_mesh(mesh, color='white', opacity=.3)
                    # p.add_mesh(scar_mesh, color='red')
                    # p.show()
                    tree = nb.KDTree(mesh.points)
                    dist, nbs = tree.query(scar_mesh.points, k=k)
                    nbs = list(set(np.array(nbs).flatten()))
                close_points[name] = nbs

    if len(auto_detected):  # not all regions were manually selected: some were present and functional in .mesh file
        if case in ['OC52']:
            tagdict = {'MV': [-8], 'LPV': [-5, -4], 'RPV': [-6, -7]}
        else:
            tagdict = {'MV': [-4], 'LPV': [-5], 'RPV': [-6]}
        meshToCsv(glob.glob(d+'*.mesh')[0], meshdir=d, verbose=False)  # make tagged vtk file

        cm = colorFromCsv(d)  # tagged (color) mesh
        for name in auto_detected:  # tags present in CARTO .mesh file
            t = tagdict[name]  # tag number
            # indices of colored mesh that have tag t: scar
            scar = [ind for ind, tag in zip(range(len(cm.points)), cm['color']) if tag == t]
            # For each mesh point, find closest colored point
            if '.vtk' in meshname:
                tree = nb.KDTree(cm.points)
                _, color_points = tree.query(mesh.points, k=1)
                # color_points is now a list of length len(mesh.points) and each value refers to an
                # index of colored mesh
                color_points = np.array(color_points).flatten()
                # colored mesh indices to vtk mesh indices
                # only consider mesh points now that are close to colored point with tag t
                scar = [i for i, e in enumerate(color_points) if e in scar]
                # Split mesh into scar region and the rest
                ind = [e for e in range(len(mesh.points)) if e not in scar]
                scar_ind = [e for e in range(len(mesh.points)) if e in scar]
                scar_ind = random.sample(scar_ind, len(scar_ind)//10)
                scar_mesh = pv.PolyData(mesh.points[scar_ind])
                cut_mesh = pv.PolyData(mesh.points[ind])  # cut away scar
                cut_mesh['Id'] = ind  # save original indexing
                tree = nb.KDTree(cut_mesh.points)
                _, nbs = tree.query(scar_mesh.points, k=k)
                nbs = list(set(np.array(nbs).flatten()))  # filter out double indices
                nbs = cut_mesh['Id'][nbs]
            elif '.csv' in meshname:
                # csv files (based on simulation data) don't have scar regions in the mesh (no activation recorded)
                # this gives a significant speedup
                # flip around the two trees: f/e colored point, find k closest mesh points
                # immediately find closest neighbors for the mesh
                tree = nb.KDTree(mesh.points)
                _, nbs = tree.query(cm.points[scar], k=k)
                nbs = list(set(np.array(nbs).flatten()))  # filter out double indices
                # p = pv.Plotter()
                # p.add_mesh(mesh, color='white', opacity=.3)
                # p.add_mesh(pv.PolyData(cm.points[scar]), color='red')
                # p.show()

            # p = pv.Plotter()
            # p.add_mesh(mesh, color='white', opacity=.4, label='Mesh')
            # p.add_mesh(pv.PolyData(mesh.points[nbs]), color=colors[2], label='Neighbors')
            # # p.add_mesh(scar_mesh, color=colors[1], label='Scar')
            # p.add_legend()
            # p.show()
            close_points[name] = nbs
    return OrderedDict(sorted(close_points.items()))


def transform(case='OC63', location='../carto_meshes', reference_case='OC45', plot=False, return_transformation=False, mesh=None,
              reference_mesh=None, reference_cores=None):
    """
    Takes input mesh, calculates the cores of the anatomical object (LPV, RPV, MV), calculates the transformation
    matrix between these three points and applies this transformation to all the mesh points.
    :param case: input mesh to be transformed. Meshname is found automatically based on case number.
    :param reference_case: name of reference mesh. Output mesh will look like this. Found automatically based on
    case number
    :param return_rot: if True, returns the Rotation object itself instead of the rotated mesh
    """

    if not mesh:
        mesh_ = pv.read(glob.glob(location + '/' + case + '/*µm.1.vtk')[0])
    else:
        mesh_ = mesh

    cores = list(getRegionCenters(case, method='median', location=location).values())
    cores -= np.mean(cores, axis=0)
    if reference_mesh is None:
        refmesh = pv.read(glob.glob(location + '/'+reference_case+'/*µm.1.vtk')[0])
    else:
        refmesh = reference_mesh
    if reference_cores is None:
        refcores = list(getRegionCenters(reference_case, method='median', location=location).values())
        refcores = np.array([p - np.mean(refcores, axis=0) for p in refcores])
    else:
        refcores = reference_cores

    # transform based on anatomical region centers
    R, _ = Rotation.match_vectors(refcores, cores, weights=[1, 2, 1])
    m = np.mean(mesh_.points, axis=0)
    if return_transformation:
        return R, m
    else:
        mesh_ = pv.PolyData(np.array([p - m for p in mesh_.points]), mesh_.cells)
        m2 = np.mean(refmesh.points, axis=0)
        refmesh.points = np.array([np.array([p - m2])[0] for p in refmesh.points])
        transf_cores = R.apply(cores)
        transf_mesh = pv.PolyData(R.apply(mesh_.points), mesh_.faces)
        if plot:
            p = pv.Plotter()
            p.add_mesh(transf_mesh, opacity=.2, color=colors[0], label='Aligned mesh: {}'.format(case),
                       show_edges=False)
            p.add_mesh(pv.PolyData(transf_cores), color=colors[0], point_size=25, label="Aligned cores")
            p.add_mesh(mesh_, opacity=.05, color=colors[1], label='Unaligned mesh: {}'.format(case), show_edges=False)
            p.add_mesh(pv.PolyData(cores), color=colors[1], point_size=25, label="Unaligned cores")
            p.add_mesh(refmesh, opacity=.05, color='black', label='Reference mesh: {}'.format(reference_case),
                       show_edges=False)
            p.add_mesh(pv.PolyData(refcores), color='black', point_size=25, label="Reference cores")
            p.background_color = "f0f0f0"
            p.add_legend(bcolor=(240/255, 240/255, 240/255), size=(.25, .25))

            # p.add_text('{} compared to {}'.format(case, reference_case))
            p.show()

        return transf_cores, transf_mesh


if __name__ == '__main__':
    transform(case='OC4', location='/media/bjorge/BACKUP/carto_meshes', plot=True)
