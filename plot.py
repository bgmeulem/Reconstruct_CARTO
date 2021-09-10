import random
import matplotlib.colors
import pyvista as pv
import scipy.stats
from create_surface import pvToPmCells
import glob
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import sys
from sklearn import neighbors as nb
import carto2csv
import add_points
sys.path.append('/'.join(os.getcwd().split('/')[:-1]))  # so that other Directories can be imported
from mesh_tools import pmToPvFaces
import create_surface
plt.style.use('fivethirtyeight')
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  # six 'fivethirtyeight' themed colors


def run(filename, title="", save=False, show_edges=True, wfcolor='white', wfopacity=.2, legend=False, half=False):
    """Plots mesh as half wireframe and half filled cells. Cell color can be set to the cell quality.
    @param filename: <str> name of the .vtk file to be plotted
    @param title: <str> title of the plot
    @param save: <bool> save the plot
    @param show_edges: <bool> show edges of mesh
    @param wfcolor: <str> color of the wireframe
    @param wfopacity: <float> Opacity of the wireframe plot (between 0-1)
    @param legend: <bool> add legend or not
    @param half: <bool> show plot as half wireframe, half filled cells"""
    if type(filename) == str:
        data = pv.read(filename)
        title = filename
    else:
        data = filename
    cells = pvToPmCells(data.cells)
    cell_center = data.points[cells].mean(1)
    cx, cy, cz = data.points.mean(axis=0)  # center of entire mesh
    quality = data.compute_cell_quality(quality_measure='volume')
    plotter = pv.Plotter(title=title)
    quality.plot(scalars="CellQuality")

    if half:
        # extract cells below the 0 xy plane
        mask = cell_center[:, 0] < cx
        mask2 = cell_center[:, 0] >= cx
        cell_ind = mask.nonzero()[0]
        subgrid = data.extract_cells([cell_ind])
        cell_ind2 = mask2.nonzero()[0]
        subgrid2 = data.extract_cells([cell_ind2])

        # advanced plotting
        plotter.add_mesh(subgrid, lighting=True, show_edges=show_edges)
        plotter.add_mesh(subgrid2, wfcolor, 'wireframe', opacity=wfopacity)
        if legend:
            plotter.add_legend([[' Added points ', wfcolor],
                                [' Tesselated Mesh ', 'black']])
    else:
        plotter.add_mesh(data, lighting=True, show_edges=show_edges)
    if save:
        plotter.show(screenshot=filename.split(".")[0] + '.png')
    else:
        plotter.show()


def plotStimBlock(location='/media/bjorge/BACKUP/carto_meshes', case='OC63', stim=4):
    meshname = glob.glob('{}/{}/*µm.1.vtk'.format(location, case))[0]
    vtkmesh = pv.read(meshname)
    s, b = pd.read_csv('{}/{}/StimBlocks/stim{}.csv'.format(location, case, stim)),\
           pd.read_csv('{}/{}/StimBlocks/block{}.csv'.format(location, case, stim))
    rpv, lpv = pd.read_csv('{}/{}/Scars/RPV.csv'.format(location, case)),\
               pd.read_csv('{}/{}/Scars/LPV.csv'.format(location, case))
    p = pv.Plotter()
    p.add_mesh(vtkmesh, color='lightgrey')
    p.add_mesh(vtkmesh.points[s.values.flatten()], color=colors[5], label='Stimulus')
    p.add_mesh(vtkmesh.points[b.values.flatten()], color='black', label='Block')
    p.add_mesh(vtkmesh.points[lpv.values.flatten()], color=colors[0], label='LPV')
    p.add_mesh(vtkmesh.points[rpv.values.flatten()], color=colors[2], label='RPV')
    p.background_color = 'f0f0f0'
    p.add_legend(bcolor=(240/255, 240/255, 240/255))
    p.show()


def writeBeforeAfter(meshdir='', refine_steps=7, boxplot=False, plot_mesh=False):
    """Writes out a file of the mesh before and after refinement"""
    meshfile = glob.glob("{}/*.mesh".format(meshdir))[0]
    print(meshfile)
    carto2csv.cartoToCsv(meshfile, meshdir)
    if not os.path.exists('{}/endo.txt'.format(meshdir)):
        add_points.writePoints(meshdir)
    endo = pd.read_csv(meshdir + 'endo.txt', sep=',', index_col=0).values
    epi = pd.read_csv(meshdir + 'epi.txt', sep=',', index_col=0).values
    triangles = pd.read_csv(meshdir + 'TrianglesSection.csv', sep=',')
    indices = triangles[['Vertex0', 'Vertex1', 'Vertex2']].values
    before_mesh = pv.PolyData(endo, pmToPvFaces(indices))
    pv.save_meshio('{}{}'.format(meshdir, 'before_mesh.vtk'), before_mesh)
    after_mesh = create_surface.run('after_mesh', meshdir, boxplot=False, plot_mesh=plot_mesh,
                                    refine_steps=10, returnmesh=True,
                                    max_edge=1200, min_edge=400)
    # pv.plot(after_mesh)
    inbetween_mesh = before_mesh + pv.PolyData(epi, pmToPvFaces(indices))
    pv.save_meshio(meshdir+'inbetween_mesh.vtk', inbetween_mesh)
    pv.save_meshio(meshdir+'after_mesh.vtk', after_mesh)
    if plot_mesh:
        p = pv.Plotter()
        p.add_mesh(before_mesh, color='f0f0f0')
        p.background_color = (255, 255, 255)
        p.update_scalar_bar_range((0, 4))
        p.show()


def plotCVVariationMethod(meshdir=''):
    meshfile = glob.glob("{}/*µm_CV0.vtk".format(meshdir))[0]
    mesh = pv.read_meshio(meshfile)
    mesh = mesh.ctp()
    df = pd.DataFrame({"index": range(mesh.n_points), "speed": mesh["speed"]})
    p = pv.Plotter()
    p.background_color = "#f0f0f0"

    points = mesh.points
    tree = nb.KDTree(points)
    i = random.randint(0, len(points))
    print(i)
    # i = 68657  # light blue
    # i = 25073  # pink
    i = 22019  # bit of both
    pt = points[i]
    dist, surfacenbs = tree.query([pt], k=10)
    # i = 18749

    dist, neighbors = tree.query([pt], k=100)
    neighborCVs = df.loc[[int(e) for e in neighbors[0]]]["speed"]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [colors[0], (240 / 255, 240 / 255, 240 / 255), colors[1]])
    sargs = dict(
        title_font_size=28,
        label_font_size=24,
        n_labels=3,
        fmt="%.1f",
        bold=True,
        font_family="arial",
        color='black'
    )
    p.add_mesh(mesh, cmap=cmap, scalars=mesh["speed"])
    p.add_mesh(mesh.points[surfacenbs[0]][5], color=colors[2], point_size=40)
    p.add_mesh(pv.PolyData(mesh.points[neighbors]), color='f0f0f0', point_size=20)

    p.show()

    freq, bins, patches = plt.hist(neighborCVs, density=True)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    # scale values to interval [0,1]
    col = bin_centers - min(bin_centers)
    col /= max(col)
    for speed, p in zip(bin_centers, patches):
        plt.setp(p, 'facecolor', cmap(speed/1.4))
        plt.setp(p, 'edgecolor', 'grey')
    mean, sigma = np.mean(neighborCVs), np.std(neighborCVs)
    plt.vlines(mesh["speed"][i], color=colors[2],ymin=0, ymax=scipy.stats.norm.pdf(mesh["speed"][i], mean, sigma),
               label='Original CV')
    new_cv = np.random.normal(mean, sigma, 1)
    plt.vlines(np.abs(new_cv), color=colors[3], ymin=0, ymax=scipy.stats.norm.pdf(new_cv, mean, sigma),
               label='New CV')
    x = np.linspace(bins[0], bins[-1], 100)
    plt.plot(x, scipy.stats.norm.pdf(x, mean, sigma), color='darkgrey')
    plt.tick_params(axis='y', left=False, labelleft=False)
    plt.xlabel('Conduction velocity (mm/ms)')
    plt.legend()
    plt.title("Conduction velocity variation")
    plt.tight_layout()
    plt.show()
    plt.savefig('../Plots/CVVariation.png', dpi=300)


if __name__ == '__main__':
    # run("../carto_meshes/OC4/1-AT1260-CLEANED.vtk")
    plotStimBlock(case='OC63', stim=2)
    # writeBeforeAfter(meshdir='/media/bjorge/BACKUP/carto_meshes/OC63/')
    # plotCVVariationMethod('/media/bjorge/BACKUP/carto_meshes/OC63')
