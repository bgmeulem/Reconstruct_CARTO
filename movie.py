from carto_mesh import *


def getSubMesh(mesh, n_, side="left"):
        x_min, x_max, y_min, y_max, z_min, z_max = mesh.bounds
        cell_centers = mesh.cell_centers().points
        if side == "left":
            mask = cell_centers[:, 0] < x_min + n_ * (x_max - x_min)
        elif side == 'right':
            mask = cell_centers[:, 0] >= x_min + n_ * (x_max - x_min)
        cell_ind = mask.nonzero()[0]
        submesh = mesh.extract_cells(cell_ind)
        return submesh, submesh.active_scalars


def plot2Meshes(mesh1: pv.PolyData, mesh2: pv.PolyData, mesh1_params, mesh2_params, n: float, name: str):
    p = pv.Plotter(off_screen=True)
    s1, scalars = getSubMesh(mesh1, n, 'left')
    s2, scalars2 = getSubMesh(mesh2, n, 'right')
    p.add_mesh(s1, **mesh1_params)
    p.remove_scalar_bar()
    p.add_mesh(s2, **mesh2_params)
    p.remove_scalar_bar()
    p.update()
    p.show(screenshot='Frames/{}_{:.2f}.png'.format(name, 1-n))


def writeFrames(mesh1, mesh2, mesh1_params, mesh2_params, name):
    for n in tqdm(np.linspace(0.01, .99, 100), desc="Writing frames"):
        plot2Meshes(mesh1, mesh2, mesh1_params, mesh2_params, n, name)


def writeFramesCartoToDoubleLayer(meshname: str = 'BlankMeshes/OC59'):
    mesh1 = CartoMesh(meshname)
    mesh2 = CartoMesh(meshname)
    mesh2.splitLayer()
    mesh1_params = {}
    mesh2_params = {}
    writeFrames(mesh1.mesh, mesh2.mesh.extract_all_edges(), mesh1_params, mesh2_params, 'cartoToDoubleLayer')


def writeFramesRefine(meshname: str = 'BlankMeshes/OC59'):
    mesh1 = CartoMesh(meshname)
    mesh1.splitLayer()
    mesh2 = CartoMesh(meshname)
    mesh2.splitLayer()
    mesh2.homogenizeMesh()
    mesh1_params = {'clim': [min(mesh1.mesh.active_scalars), 0], 'below_color': 'blue', 'above_color': 'red'}
    mesh2_params = {'color': 'white', 'opacity': .3}
    writeFrames(mesh1.mesh.extract_all_edges(), mesh2.mesh.extract_all_edges(), mesh1_params, mesh2_params, 'refine')


def writeFramesTetrahedralise(meshname: str = 'BlankMeshes/OC59'):
    mesh1 = CartoMesh(meshname)
    mesh1.splitLayer()
    mesh1.homogenizeMesh()
    mesh2 = pv.read('BlankMeshes/OC59/OC59_MV_only_600-1000µm.1.vtk')
    mesh1_params = {'color': 'white', 'opacity': .3}
    mesh2_params = {}
    writeFrames(mesh1.mesh.extract_all_edges(), mesh2, mesh1_params, mesh2_params,
                'tetrahedralise')


def writeFramesConductionVelocity():
    mesh1 = pv.read('BlankMeshes/OC59/OC59_MV_only_600-1000µm.1.vtk')
    mesh1_params = {}
    mesh2 = pv.read('BlankMeshes/OC59/OC59_MV_only_600-1000µm_CV0.vtk')
    mesh2_params = {}
    writeFrames(mesh1, mesh2, mesh1_params, mesh2_params,
                'cv')


if __name__ == '__main__':
    # writeFramesCartoToDoubleLayer()
    # writeFramesRefine()
    # writeFramesTetrahedralise()
    writeFramesConductionVelocity()
