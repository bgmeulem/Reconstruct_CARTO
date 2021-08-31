import pandas as pd
import numpy as np
import pyvista as pv
from tqdm import tqdm
from mesh_tools import pvToPmFaces


def normalize(vector):
    return [e / np.sqrt(np.dot(vector, vector)) for e in vector]


def writePoints(meshdir="", from_stl=False):
    if from_stl:
        mesh = pv.read("test.stl")
        mesh.compute_normals(inplace=True)
        vert = pd.DataFrame(mesh.points, columns=[['X', 'Y', 'Z']])
        tri = pd.DataFrame(pvToPmFaces(mesh.faces), columns=[['Vertex0', 'Vertex1', 'Vertex2']])
        tri[["NormalX", "NormalY", "NormalZ"]] = pd.DataFrame(mesh["Normals"], index=tri.index)
        tri["GroupID"] = mesh.n_faces*[0]  # ignore GroupID
        tri.to_csv("TrianglesSection.csv")
        vert.to_csv("VerticesSection.csv")
    else:
        vert = pd.read_csv(meshdir+"VerticesSection.csv", sep=',')
        tri = pd.read_csv(meshdir+'TrianglesSection.csv', sep=',')
    endo = open(meshdir+'endo.txt', 'w+')
    endo.write(meshdir+'Index,X,Y,Z\n')
    epi = open(meshdir+'epi.txt', 'w+')
    epi.write(meshdir+"Index,X,Y,Z\n")
    mid = open(meshdir+'mid.txt', 'w+')
    mid.write(meshdir+"Index,X,Y,Z\n")
    # allpoints = open('allpoints.txt', 'w+')
    # allpoints.write('Index,X,Y,Z\n')
    tqbar = tqdm(range(vert.shape[0]), desc='        Adding points')
    for index, point in vert.iterrows():
        tqbar.update(1)
        co = point[['X', 'Y', 'Z']]

        # Select all facets that uses the vertex (usually each vertex appears in  ~5 triangles)
        facets = tri.loc[(tri[['Vertex0']].values == index) | (tri[['Vertex1']].values == index) |
                         (tri[['Vertex2']].values == index)]
        av_normal = [0, 0, 0]

        # loop over all facets and calculate average normal
        for row in facets[['NormalX', 'NormalY', 'NormalZ', 'GroupID']].iterrows():
            for i in range(len(av_normal)):  # from 0 to 2
                av_normal[i] += row[1][i]
        av_normal = normalize(av_normal)
        newpoint = [e + 0.4 * n for e, n in zip(co, av_normal)]
        newpoint2 = [e - 0.1 * n for e, n in zip(co, av_normal)]

        endo.write(str(index))
        epi.write(str(index))
        mid.write(str(index))
        for co1, co2, co3 in zip(co, newpoint, newpoint2):
            # Convert from mm to µm
            mid.write(',' + str(co1 * 1000.))
            epi.write(',' + str(co2 * 1000.))
            endo.write(',' + str(co3 * 1000.))

        epi.write('\n')
        mid.write('\n')
        endo.write('\n')

    endo.close()
    epi.close()
    mid.close()
    tqbar.close()
    print("")


if __name__ == '__main__':
    writePoints()
