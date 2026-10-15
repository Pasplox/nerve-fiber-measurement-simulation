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
loc_x = 50 # cm
loc_y = 0.3 # cm

# extracellular properties
sigma_o = 5 # mS/cm | extracellular conductivity

# axonal properties
a = 300 # um | axon radius
v_c = 22.5 # mV/s | conduction velocity