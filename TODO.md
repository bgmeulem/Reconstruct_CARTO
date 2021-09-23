# Requirements
- check if it works with PyMesh and Triangle libraries
- Move all dependencies to one requirements.txt file
- make python script setup.py to pip install requirements.txt and install TetGen, Triangle and PyMesh

# Cleanup
- Big cleanup: re-order everything to match perhaps a class? Think about better structure
  - One class per step, each time inheriting from the previous one? Maybe not very memory-efficient
- Remove anything unnecessary

# Implementation
- Add DGM command to calculate velocities

# Extra tips
- Use dictionaries to avoid overuse of if-else statements:
  - add cases to consider as dictionary keys and resulting behaviour (function, object, value...) as dict value
- optional parser arguments as addvalue instead of const-default (which looks like they all need to have a default value)
- docstring in sfinx style instead of doxygen style. Sfinx style looks something like (look this up):
  - Args:
    - args
    - kwargs
  - Returns
    - (type): returns type (type) containing:
      - value1
      - value2