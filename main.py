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
xp = 50 # cm | electrode x position (from lab handout)
Ro = 0.3 # cm | electrode y position / radial distance from fiber (from lab handout)

# extracellular properties
sigma_o = 5 # mS/cm | extracellular conductivity (from lab handout)

# intracellular properties
resistivity_i = 125 # ohm-cm | intracellular resistivity (from Cable Eqn Params)
sigma_i_Scm = 1 / resistivity_i # S/cm | intracellular conductivity
sigma_i = sigma_i_Scm * 1000 # mS/cm | intracellular conductivity

# axonal properties
a_um = 300 # um   | axon radius (from lab handout)
a = a_um * 1e-4 # cm   | axon radius
v_c_ms = 22.5 # m/s  | conduction velocity (from lab handout)
v_c = v_c_ms * 100 / 1000 # cm/ms | conduction velocity

# membrane properties
Rm = 3000 # ohm-cm^2 | membrane resistance (from Cable Eqn Params)
Cm = 1 # uF/cm^2 | membrane capacitance (from Cable Eqn Params)

# fiber geometry
L = 1000 # cm   | fiber length
dx = 100 * (2 * a) # cm   | node spacing = 100 * axonal diameter (approx from class)
xj = np.arange(0, L, dx) # cm   | node positions along fiber
N = len(xj) # number of nodes

# time parameters
dt = 0.01 # ms   | time step
T = 50 # ms   | total simulation time
t = np.arange(0, T, dt) # ms   | time array

# action potential parameters
V_rest = -70 # mV   | resting membrane potential
V_AP = 100 # mV   | AP amplitude above rest

# gaussian AP model parameters
AP_duration = 1 # ms   | approximate AP duration at a point
AP_FWHM = AP_duration * v_c # cm   | spatial FWHM (duration × conduction velocity)
AP_sigma = AP_FWHM / 2.36 # cm   | gaussian spatial sigma

### calculated parameters
coeff = -(a**2 * sigma_i) / (4 * sigma_o) # cm^2 | dipole coefficient (-a^2 * sigma_i / 4 * sigma_o)
# note: we only have to calculate these once since the electrode position is constant
rj = np.sqrt((xj - xp)**2 + Ro**2)        # cm   | distance from each node to electrode
d = (xp - xj) / rj**3                     # cm^-2 | geometric distance term (N x 1)
# to calculate B matrix, format matrix to calculate central difference derivative of V_m
# first and last rows remain zero (boundary condition)
B = np.zeros((N, N))
diag_idx = np.arange(1, N - 1)
B[diag_idx, diag_idx - 1] = -1 / (2 * dx) # subdiagonal
B[diag_idx, diag_idx + 1] =  1 / (2 * dx) # superdiagonal