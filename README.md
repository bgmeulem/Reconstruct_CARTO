# Reconstruct_CARTO
Reconstructs a simulatable 3D mesh from clinical CARTO mapping data, providing control over mesh resolution, 
conduction velocity distributions and non-conductive regions.

<p align="center">
<img src=https://media0.giphy.com/media/bZM2OaOQb4HCVymzna/giphy.gif?cid=790b76112533f1a2f99b476d6833aa55d4b4c8ef9e3227b2&rid=giphy.gif&ct=g />
</p>

# Requirements:
Some requirements are not automatically installable. These requirements are:
- [tetgen 1.6.0](http://www.wias-berlin.de/software/index.jsp?id=TetGen&lang=1#Download) ([Download here](http://www.wias-berlin.de/software/tetgen/download2.jsp))
  - The ```tetgen``` command should be aliased to the tetgen executable, located at
```<tetgen_installation_folder>/tetgen1.6.0/build/tetgen```
- [PyMesh 0.3](https://pymesh.readthedocs.io/en/latest/index.html)  ([Download here](https://pymesh.readthedocs.io/en/latest/installation.html))
  - Do **NOT** install PyMesh with pip. This is another module with the same name. Download, build and install according to the setup as provided by the PyMesh docs.
  - Installing the third-party dependencies can be done as specified in the docs, but also by running ```./build.py all``` in the directory ```third_party```.


To install the remaining requirements, run
```
pip install -r requirements.txt
```
This will install the following modules, or update them if they are not already installed:
- [PyVista](https://docs.pyvista.org/getting-started/index.html)
- tqdm
- matplotlib
- numpy
- sklearn
- pandas


# Usage
**Reconstructing** involves taking an input ```.mesh``` file and making it into a tetrahedron mesh with user-defined 
resolution. Conduction velocities can already be interpolated on this mesh if a .csv file containing point coordinates 
and speed values is passed as an optional argument.
```
python reconstruct.py <name> --<speed_file>
```
**Interpolating conduction velocities** can also be done after reconstruction, e.g. in case you want to select
non-conductive regions on the mesh *after* the reconstruction. To this end, a .csv file with point coordinates and
conduction velocity values need to be passed as an argument to ``apply_cv.py``.
```
python apply_cv.py <name> --<write_adjust> --<region_dir> --<speed_file> --<ncv> --<speed_col> --<writeVTK>
```

Alternatively, if you want full control over the mesh reconstruction and its subfacets, you can also import the
CartoMesh class and its dependencies.
```
>>> from carto_mesh import *
>>> m = CartoMesh('filename.mesh')
>>> m.reconstruct()
```

# Documentation
Documentation can be found [here](https://bgmeulem.github.io/Reconstruct_CARTO/html/index.html)
