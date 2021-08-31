import pyvista as pv
import glob
import reconstruct_mesh
import pandas as pd
import numpy as np
import create_surface
import pymesh as pm
from pymeshfix import MeshFix


def printDistances(segments, mesh):
    for i in range(len(segments) // 2):
        ids = segments[2 * i] + segments[2 * i + 1]
        unique_ids = []
        for e in ids:
            if ids.count(e) > 1:
                commonpoint = mesh.points[e]
                unique_ids = [id for id in ids if id != e]
        co1, co2 = mesh.points[unique_ids[0]], mesh.points[unique_ids[1]]
        print("Pair ", i)
        print(create_surface.dist(commonpoint, co2))
        print(create_surface.dist(commonpoint, co1))


def run(meshname):
    """This is a quick n dirty aiding tool for when TetGen fails. Printing out the failed edges
    in mesh_tools.makeNonColinear() allows for a direct copy paste from the terminal in here.
    This function plots the mesh after the refinement procedure and colors all the points of the
    failed edges in red. Dirty and hardcoded, but better than nothing."""
    segments = [[[352312, 352284], [352312, 428493]], [[22475, 22484], [22475, 22494]],
                [[421930, 220231], [421930, 421931]], [[221439, 221440], [221439, 428339]],
                [[206200, 206210], [206200, 206211]], [[205156, 205157], [205156, 428222]],
                [[402643, 423333], [402643, 423334]], [[428168, 190625], [428168, 428169]],
                [[113442, 113434], [113442, 204981]], [[420500, 235296], [420500, 428748]],
                [[78883, 78893], [78883, 78898]], [[419442, 419439], [419442, 428669]],
                [[424519, 291217], [424519, 428412]], [[140132, 140141], [140132, 140148]],
                [[316816, 388877], [316816, 388876]], [[388899, 352284], [388899, 428494]],
                [[248998, 249001], [248998, 249004]], [[248914, 423618], [248914, 423621]],
                [[204977, 113418], [204977, 427984]], [[407377, 408776], [407377, 428603]],
                [[172334, 136998], [172334, 172333]], [[127613, 50510], [127613, 167315]],
                [[100899, 100891], [100899, 427941]], [[423335, 423329], [423335, 428707]],
                [[231338, 305706], [231338, 305709]], [[305679, 305677], [305679, 305682]],
                [[946, 86299], [946, 86300]], [[100320, 206948], [100320, 427938]],
                [[154035, 185094], [154035, 185095]], [[423333, 423331], [423333, 423332]],
                [[352284, 352283], [352284, 428496]], [[78744, 78740], [78744, 78747]],
                [[204983, 113419], [204983, 427982]], [[205469, 195951], [205469, 205474]],
                [[218959, 428328], [218959, 428330]], [[113444, 204983], [113444, 427979]],
                [[255335, 420486], [255335, 420485]], [[426007, 100890], [426007, 426008]],
                [[242701, 242692], [242701, 252125]], [[5497, 427846], [5497, 427847]],
                [[297545, 423763], [297545, 423762]], [[325934, 325928], [325934, 428459]],
                [[421624, 421622], [421624, 428682]], [[173651, 114721], [173651, 173652]],
                [[380351, 380354], [380351, 380355]], [[305016, 349854], [305016, 428490]],
                [[260437, 260433], [260437, 428409]], [[173669, 95099], [173669, 428127]],
                [[380396, 380417], [380396, 380418]], [[222328, 222331], [222328, 222329]],
                [[185053, 185050], [185053, 185052]], [[78740, 78742], [78740, 78743]],
                [[140438, 428027], [140438, 428026]], [[320936, 407454], [320936, 428597]],
                [[412695, 401732], [412695, 412696]], [[291268, 291266], [291268, 428421]],
                [[168333, 168344], [168333, 428104]], [[218895, 218897], [218895, 218899]],
                [[298530, 298555], [298530, 428431]], [[421614, 421606], [421614, 428678]],
                [[26882, 26929], [26882, 427877]], [[230985, 259395], [230985, 259396]],
                [[331695, 419395], [331695, 419396]], [[10648, 110252], [10648, 427973]],
                [[374297, 260259], [374297, 374301]], [[41937, 41934], [41937, 427897]],
                [[179030, 179012], [179030, 179031]], [[422837, 422844], [422837, 422855]],
                [[48912, 48907], [48912, 427909]], [[424592, 257724], [424592, 428714]],
                [[421629, 421628], [421629, 421631]], [[360773, 318013], [360773, 426993]],
                [[407360, 423387], [407360, 407356]], [[318418, 428456], [318418, 428750]],
                [[422838, 422844], [422838, 422842]], [[205992, 205970], [205992, 428243]],
                [[206998, 428045], [206998, 428752]], [[327129, 327151], [327129, 327136]],
                [[221398, 221400], [221398, 221401]], [[206226, 206227], [206226, 428270]],
                [[423228, 423235], [423228, 427730]], [[407454, 320939], [407454, 428598]],
                [[418590, 418591], [418590, 418593]], [[191886, 190550], [191886, 428189]],
                [[383813, 257227], [383813, 383848]]]  # from tetgen output
    if "/" in meshname:
        meshdir = "/".join(meshname.split("/")[:-1]) + "/"
        print("meshdir: ", meshdir)
    mesh = reconstruct_mesh.run(meshname=meshname, return_surfacemesh=True, max_edge=500, min_edge=200, refine_steps=12,
                                keep_intmed=True, skip_reading=True)

    points = mesh.points
    triangles = pd.read_csv(meshdir + "TrianglesSection.csv", sep=',')[['Vertex0', 'Vertex1', 'Vertex2']].to_numpy()

    plotter = pv.Plotter(title="{} pairs of intersecting segments".format(len(segments) / 2))
    plotter.add_mesh(mesh, show_edges=True, opacity=.3)
    for edgepair in segments:
        edge1, edge2 = edgepair[0], edgepair[1]
        mesh = pv.PolyData(points[edge1+edge2])
        plotter.add_mesh(mesh, color="red", point_size=20)
    plotter.show()


if __name__ == '__main__':
    dirname = 'carto_meshes/OC49/'
    meshname = glob.glob('../' + dirname + '*.mesh')[0]
    run(meshname)
