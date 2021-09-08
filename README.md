# Reconstruct_CARTO
Reconstructs a simulatable 3D mesh from CARTO data

reconstruct_mesh.py combines the following files to reconstruct a simulatable mesh from carto data:
1. carto2csv.py
2. add_points.py
3. create_surface.py
4. tet.py
5. mesh_tools.py
6. apply_cv.py

If you want to apply conduction velocities, you should have a .csv file in the same
directory as the mesh called "speed.csv" with the calculated speeds and column names "x,y,z,speed".

TODO:
- not all intermediate files get deleted
- implementation of reading in .stl meshes
- implementation of reading in .mat meshes
- Automate calculation of conduction velocities with DGM
- Automatically generate "speed.csv" with DGM

