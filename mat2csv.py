import numpy as np
import scipy.io
import pyvista as pv
from mesh_tools import pvToPmFaces, pmToPvFaces
from tqdm import tqdm
import time
import sklearn.neighbors
import itertools


def read_matlab_file(filename):
    data = scipy.io.loadmat('{0}'.format(filename))['data']
    # maps = data[0][0]['maps']
    # activation = maps[0][0]['activation']
    # bipolar_ms = activation[0][0]['bipolar_ms']
    # unipolar_ms = activation[0][0]['unipolar_ms']
    # voltage = maps[0][0]['voltage']
    # bipolar_mV = voltage[0][0]['bipolar_mV']
    # unipolar_mV = voltage[0][0]['unipolar_mV']
    anatomy = data[0][0]['anatomy'][0][0]
    vertices = anatomy['vertices_mm']
    # faces = anatomy['faces']  # incorrect

    # xyz, bip = DGMdefinitions.matlab_to_coord(0.05, bipolar_ms, unipolar_ms, bipolar_mV, unipolar_mV, vertices, maps,
    #                                          period=200)
    # self.length = len(xyz)
    # return xyz, bip, vertices, maps
    return vertices


def generateNormals(vertices, k=5):
    """Generates normals of vertices based on k closest points"""
    c = np.mean(vertices, axis=0)  # center of mesh
    tree = sklearn.neighbors.KDTree(vertices, leaf_size=2 * k)
    dist, ind = tree.query(vertices, k=k)
    normals = []
    for i, point in enumerate(vertices):
        radial = [p - co for p, co in zip(point, c)]
        nbs = ind[i]  # indices of k neighbours
        totnormal = [0, 0, 0]
        for nb1, nb2 in list(itertools.combinations(nbs, 2)):
            normal = np.cross(vertices[nb1], vertices[nb2])
            if np.dot(normal, radial) < np.pi / 2:
                normal = [-e for e in normal]
            totnormal = [sum(e) for e in zip(totnormal, normal)]
        totnormal = [e / np.sqrt(np.dot(totnormal, totnormal)) for e in totnormal]
        normals.append(totnormal)
    return normals


def writePly(points, normals, name):
    with open(name.split('.')[0] + '.ply', 'w+') as of:
        # Header
        of.write("ply\nformat ascii 1.0\nelement vertex {}\n"
                 "property float x\nproperty float y\nproperty float z\n"
                 "property float nx\nproperty float ny\nproperty float nz\n"
                 "element face 0\nproperty list uchar int vertex_index\n"
                 "end_header".format(len(points)))
        for point, normal in zip(points, normals):
            of.write('\n')
            for c in point:
                of.write(str(c) + ' ')
            for n in normal[:2]:
                of.write(str(n) + ' ')
            of.write(str(normal[2]))
    of.close()


def writeCsv(points, faces, meshdir):
    tq = tqdm(range(len(points)), desc='        Writing points')
    with open(meshdir+'VerticesSection.csv', 'w+') as of:
        of.write("Ind,X,Y,Z\n")

        for i in tq:
            point = points[i]
            time.sleep(1e-6)
            of.write(str(i))
            for comp in point:
                of.write("," + str(comp))
            of.write("\n")
            tq.update(1)
    tq.close()
    of.close()
    time.sleep(1e-6)

    cx, cy, cz = c = np.mean(points, axis=0)  # center of mesh
    tq = tqdm(range(len(faces)), desc='        Writing faces', )
    with open(meshdir+"TrianglesSection.csv", "w+") as of:
        of.write("Index,Vertex0,Vertex1,Vertex2,NormalX,NormalY,NormalZ,GroupID\n")

        for i in tq:
            time.sleep(1e-6)
            face = faces[i]
            # calculate normal
            p0, p1, p2 = points[face]
            ct = np.mean([p0, p1, p2], axis=0)  # center of triangle
            x0, y0, z0 = p0
            x1, y1, z1 = p1
            x2, y2, z2 = p2
            # plane orientation defining vectors:
            ux, uy, uz = [x1 - x0, y1 - y0, z1 - z0]  # p1 -> p0 = vec(u)
            vx, vy, vz = [x2 - x0, y2 - y0, z2 - z0]  # p2 -> p0 = vec(v)

            normal = [uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx]  # unnormalized, happens later

            # check if normal is correct direction
            radial = [a - b for a, b in zip(ct, c)]  # from center of mesh to center point of triangle
            if np.dot(normal, radial) < np.pi / 2:  # normal points inwards
                normal = [-e for e in normal]  # now it points outwards

            of.write(str(i) + ",")  # face index
            for vert in face:  # face vertex indices
                of.write(str(vert) + ",")
            for comp in normal:
                of.write(str(comp) + ",")
            of.write(str(0) + '\n')  # GroupID (unused)
            tq.update(1)
    tq.close()
    of.close()
    time.sleep(1e-6)


def mat2ply(name):
    vertices = read_matlab_file(name)

    print("        Generating normals")
    normals = generateNormals(vertices)

    print("        Writing .ply file")
    meshname = name.split('.')[0]
    writePly(vertices, normals, meshname)


def ply2csv(infile='mat_meshes/pt_2/3_LA_at2.mat'):

    meshdir = "/".join(infile.split("/")[:-2]) + "/"

    surface = pv.PolyData("mat_meshes/pt_1_test_new/meshlab.ply")
    writeCsv(surface.points, pvToPmFaces(surface.faces), meshdir)


if __name__ == '__main__':
    ply2csv(infile='mat_meshes/pt_1_test_new/1_LA_AT.mat', plot_fix=True)
