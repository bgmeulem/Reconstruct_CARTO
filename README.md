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

# Dependencies:
- tetgen 1.6.0 (version is important as output gets caught and interpreted)
- PyVista 0.27.3
- [PyMesh 0.3](https://pymesh.readthedocs.io/en/latest/installation.html)
# Python module dependencies
-tqdm
- time
- glob
- matplotlib.pyplot
- collections.OrderedDict
- configparser
- os
- numpy
- sklearn.neighbors
- pandas
- subprocess
- csv
- glob
- sys
- io
# TODO:
- Automate calculation of conduction velocities with DGM
- Automate installing dependencies
- Make two seperate runfiles: one for reconstruction and one for applying ad-hoc CV

