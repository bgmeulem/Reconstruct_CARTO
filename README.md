# Reconstruct_CARTO
Reconstructs a simulatable 3D mesh from clinical CARTO mapping data

<p align="center">
<img src=https://media0.giphy.com/media/bZM2OaOQb4HCVymzna/giphy.gif?cid=790b76112533f1a2f99b476d6833aa55d4b4c8ef9e3227b2&rid=giphy.gif&ct=g />
</p>


# Usage
```
>>> from carto_mesh import *
>>> m = CartoMesh('filename.mesh')
>>> m.reconstruct()
```

If you want to apply conduction velocities, you should have a .csv file in the same
directory as the mesh called "speed.csv" with the calculated speeds and column names "x,y,z,speed".

# Requirements:
Some requirements are not automatically installable. These requirements are:
- [tetgen 1.6.0](http://www.wias-berlin.de/software/index.jsp?id=TetGen&lang=1#Download) ([Download here](http://www.wias-berlin.de/software/tetgen/download2.jsp))
  - The ```tetgen``` command should be aliased to the tetgen executable, located at
```<tetgen_installation_folder>/tetgen1.6.0/build/tetgen```
- [PyMesh 0.3](https://pymesh.readthedocs.io/en/latest/index.html)  ([Download here](https://pymesh.readthedocs.io/en/latest/installation.html))
  - Do **NOT** install PyMesh with pip. This is another module with the same name. Download, build and install according to the setup as provided by the PyMesh docs.
  - Installing the third-party dependencies can be done as specified in the docs, but also by running ```./build.py all``` in the directory ```third_party```.


To install the remaining requirements (see below), simply run
```
pip install -r requirements.txt
```
This will install (or update) the following modules if they are not already installed:
- [PyVista](https://docs.pyvista.org/getting-started/index.html)
- tqdm
- matplotlib
- numpy
- sklearn
- pandas

# TODO:
- Automate calculation of conduction velocities with DGM
- Make two seperate runfiles: one for reconstruction and one for applying ad-hoc CV

