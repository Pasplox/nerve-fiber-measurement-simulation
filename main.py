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
xp = 50   # cm | electrode x position (from lab handout)
Ro = 0.3  # cm | electrode y position / radial distance from fiber (from lab handout)

# extracellular properties
sigma_o = 5  # mS/cm | extracellular conductivity (from lab handout)
resistivity_i = 125 # ohm-cm | intracellular resistivity (from Cable Eqn Params)
sigma_i = 1 / resistivity_i  # S/cm | intracellular conductivity

# axonal properties
a_um = 300   # um   | axon radius (from lab handout)
a = a_um * 1e-4  # cm   | axon radius 
v_c_ms = 22.5  # m/s  | conduction velocity (from lab handout)
v_c = v_c_ms * 100 / 1000  # cm/ms | conduction velocity

# membrane properties
Rm = 3000 # ohm-cm^2 | membrane resistance (from Cable Eqn Params)
Cm = 1    # uF/cm^2 | membrane capacitance (from Cable Eqn Params)
