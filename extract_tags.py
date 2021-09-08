import glob
import pyvista as pv
from mesh_tools import pmToPvFaces, colorFromCsv
import carto2csv
import csv


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


def run(meshdir="", writeVtk=False, write_scar=True):
    """Writes out .csv file with scar coördinates. Can also write out this surface as .vtk file"""
    # get mesh from GroupID's
    mesh = colorFromCsv(meshdir)
    if write_scar:
        # write tagged LPV, RPV, MV to noncond.csv and myocardium to myo.csv
        writePerTag(mesh, meshdir)
    if writeVtk:
        pv.save_meshio('colors.vtk', mesh)
    print("Tagged .vtk file written")


if __name__ == '__main__':
    # create csv files:
    print("Making .csv files\n")
    mn = glob.glob("*.mesh")[0]
    carto2csv.cartoToCsv(mn, "")
    run("", writeVtk=True)
