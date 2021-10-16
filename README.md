# Reconstruct_CARTO
Reconstructs a simulatable 3D mesh from CARTO data

<p align="center">
<img src=![hippo](https://media0.giphy.com/media/bZM2OaOQb4HCVymzna/giphy.gif?cid=790b76112533f1a2f99b476d6833aa55d4b4c8ef9e3227b2&rid=giphy.gif&ct=g />)
</p>


# Usage
```
from carto_mesh import *
m = CartoMesh('filename.mesh')
m.reconstruct()
```

If you want to apply conduction velocities, you should have a .csv file in the same
directory as the mesh called "speed.csv" with the calculated speeds and column names "x,y,z,speed".

# Dependencies:
- tetgen 1.6.0 (version is important as output gets caught and interpreted)
- If finished: DGM for calculating CVs (currently not used)
- PyVista 0.27.3
- PyMesh 0.3
# Python module dependencies
- sklearn.neighbors
- subprocess
- pandas
- glob
- os
- tqdm
- csv
- random
- sys
- scipy.spatial.transform.Rotation
- scipy.stats
- collections.OrderedDict
- numpy
- io
- matplotlib.pyplot
- argparse
- logging for DGMdefinitions (can be made optional)

# TODO:
- not all intermediate files get deleted
- implementation of reading in .stl meshes
- implementation of reading in .mat meshes
- Automate calculation of conduction velocities with DGM
- Automatically generate "speed.csv" with DGM

