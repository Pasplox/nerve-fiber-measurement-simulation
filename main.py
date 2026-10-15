"""
Lab 4: Potential due to Current Sources
BME 301 
Zach Gor
"""

### import libraries
import numpy as np
import matplotlib.pyplot as plt

### define parameters
# electrode location
xp = 50   # cm | electrode x position
Ro = 0.3  # cm | electrode y position (radial distance from fiber)

# extracellular properties
sigma_o = 5  # mS/cm | extracellular conductivity

# axonal properties
a_um = 300   # um   | axon radius
a = a_um * 1e-4  # cm   | axon radius 
v_c_ms = 22.5  # m/s  | conduction velocity
v_c = v_c_ms * 100 / 1000  # cm/ms | conduction velocity