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

# intracellular properties
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

# fiber geometry
L = 1000                      # cm   | fiber length
dx = 100 * (2 * a)           # cm   | node spacing = 100 * axonal diameter (approx from class)
xj = np.arange(0, L, dx)     # cm   | node positions along fiber
N = len(xj)                  # number of nodes

# time parameters
dt = 0.01                    # ms   | time step
T = 50                       # ms   | total simulation time
t = np.arange(0, T, dt)     # ms   | time array

