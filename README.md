# Reconstruct_CARTO
Reconstructs a simulatable 3D mesh from CARTO data

# Usage
python reconstruct.py --filename=location/of/meshfile

reconstruct_mesh.py combines the following files to reconstruct a simulatable mesh from carto data:
1. carto2csv.py
2. add_points.py
3. create_surface.py
4. tet.py
5. mesh_tools.py
6. apply_cv.py

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

