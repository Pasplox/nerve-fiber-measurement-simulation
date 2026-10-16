"""
Lab 4: Potential due to Current Sources
BME 301 
Zach Gor
"""

### import libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import t as t_dist
from sklearn.metrics import r2_score

### ---------------------------

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
L = 400 # cm   | fiber length
dx = 1 # cm   | integration node spacing
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

### ---------------------------

### calculated parameters
coeff = -(a**2 * sigma_i) / (4 * sigma_o) # cm^2 | dipole coefficient (-a^2 * sigma_i / 4 * sigma_o)
# to calculate B matrix, format matrix to calculate central difference derivative of V_m
# first and last rows remain zero (boundary condition)
B = np.zeros((N, N))
diag_idx = np.arange(1, N - 1)
B[diag_idx, diag_idx - 1] = -1 / (2 * dx) # subdiagonal
B[diag_idx, diag_idx + 1] =  1 / (2 * dx) # superdiagonal

### ---------------------------

def simulate(xp, Ro):
    # note: we only have to calculate these once since the electrode position is constant
    rj = np.sqrt((xj - xp)**2 + Ro**2)        # cm   | distance from each node to electrode
    d = (xp - xj) / rj**3                     # cm^-2 | geometric distance term (N x 1)
    # output extracellular potential
    # note: phi_o only varies with time, since we're measuring the extracellular potential at a single electrode location
    phi_o = np.zeros(len(t)) # mV | extracellular potential at electrode over time
    for i, ti in enumerate(t):
        # central position of the AP at current time ti
        x_AP = v_c * ti
        # gaussian AP model (Vm), calculated as a column vector with each entry corresponding to each node position xj
        # recalculated at each time step (due to AP propagation)
        Vm = V_rest + V_AP * np.exp(-((xj - x_AP)**2) / (2 * AP_sigma**2))
        # calculate g, the spatial derivative of Vm using the central difference method (B matrix)
        # @ = matrix multiplication operator
        g = B @ Vm
        # calculate phi_o at the current time step using the trapezoidal rule for numerical integration
        # syntax is equivalent to manually calculating and summing the area of trapezoids between each pair of points in the array
        phi_o[i] = coeff * np.trapz(d * g, dx=dx)
    return phi_o

### ---------------------------

### q1 plots
phi_o = simulate(xp, Ro)
# extracellular potential at electrode vs time
fig, ax = plt.subplots()
ax.plot(t, phi_o * 1000)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('$\\phi_o$ (µV)')
ax.set_title(f'Extracellular Potential at Electrode (x={xp} cm, y={Ro} cm) vs Time')
plt.tight_layout()
plt.savefig('q1_extracellular_potential_vs_time.png', dpi=300, bbox_inches='tight')
plt.close()

# Vm vs position at t=15ms
t_plot = 15 # ms
x_AP_plot = v_c * t_plot
Vm_plot = V_rest + V_AP * np.exp(-((xj - x_AP_plot)**2) / (2 * AP_sigma**2))
fig, ax = plt.subplots()
ax.plot(xj, Vm_plot)
ax.axvline(x=x_AP_plot, color='r', linestyle='--', label=f'AP center (x={x_AP_plot:.1f} cm)', zorder=0)
ax.legend()
ax.set_xlabel('Position (cm)')
ax.set_ylabel('$V_m$ (mV)')
ax.set_title(f'Transmembrane Voltage vs Position at t={t_plot} ms')
plt.tight_layout()
plt.savefig('q1_Vm_vs_position_t15ms.png', dpi=300, bbox_inches='tight')
plt.close()

# g vs position at t=10ms
t_plot_g = 10 # ms
x_AP_plot_g = v_c * t_plot_g
Vm_plot_g = V_rest + V_AP * np.exp(-((xj - x_AP_plot_g)**2) / (2 * AP_sigma**2))
g_plot = B @ Vm_plot_g
fig, ax = plt.subplots()
ax.plot(xj, g_plot)
ax.axvline(x=x_AP_plot_g, color='r', linestyle='--', label=f'AP center (x={x_AP_plot_g:.1f} cm)', zorder=0)
ax.legend()
ax.set_xlabel('Position (cm)')
ax.set_ylabel('$g$ = $\\partial V_m / \\partial x$ (mV/cm)')
ax.set_title(f'Spatial Derivative of $V_m$ vs Position at t={t_plot_g} ms')
plt.tight_layout()
plt.savefig('q1_g_vs_position_t10ms.png', dpi=300, bbox_inches='tight')
plt.close()

# Vm vs time at x=2cm
x_obs = 2 # cm
Vm_at_x = V_rest + V_AP * np.exp(-((x_obs - v_c * t)**2) / (2 * AP_sigma**2))
fig, ax = plt.subplots()
ax.plot(t, Vm_at_x)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('$V_m$ (mV)')
ax.set_title(f'Transmembrane Voltage vs Time at x={x_obs} cm')
plt.tight_layout()
plt.savefig('q1_Vm_vs_time_x2cm.png', dpi=300, bbox_inches='tight')
plt.close()

### ---------------------------

### q2 plots
phi_o = simulate(60, Ro)
# extracellular potential at electrode vs time
fig, ax = plt.subplots()
ax.plot(t, phi_o * 1000)
ax.set_xlabel('Time (ms)')
ax.set_ylabel('$\\phi_o$ (µV)')
ax.set_title(f'Extracellular Potential at Electrode (x={xp} cm, y={Ro} cm) vs Time')
plt.tight_layout()
plt.savefig('q2_extracellular_potential_vs_time.png', dpi=300, bbox_inches='tight')
plt.close()

### ---------------------------

### q3 plots
# determine extracellular potential (pp) at varying electrode y positions
y_values = np.arange(0.1, 5.1, 0.5)
pp_amps = np.zeros(len(y_values))  # µV
for i, Ro in enumerate(y_values):
    phi_o = simulate(50, Ro)
    pp_amps[i] = (np.max(phi_o) - np.min(phi_o)) * 1000  # µV | peak-to-peak amplitude

# save data to CSV
np.savetxt('q3_pp_amp_vs_y.csv', np.column_stack([y_values, np.round(pp_amps, 4)]),
           delimiter=',', header='y (cm),Peak-to-Peak (µV)', comments='', fmt='%.4f')

# curve fit models
# using proportionality constants A and B for monopole and dipole models, respectively
def monopole(y, A):
    return A / y

def dipole(y, B):
    return B / y**2

# fit the models to the data
popt_mono, pcov_mono = curve_fit(monopole, y_values, pp_amps)
popt_dip, pcov_dip  = curve_fit(dipole, y_values, pp_amps)

y_fit = np.linspace(y_values[0], y_values[-1], 500)

# plot peak-to-peak amplitude vs electrode-to-fiber distance with curve fits
fig, ax = plt.subplots()
ax.scatter(y_values, pp_amps, label='Data', zorder=3)
ax.plot(y_fit, monopole(y_fit, *popt_mono), label=f'Monopole fit ($A/y$, A={popt_mono[0]:.2f})', color='orange')
ax.plot(y_fit, dipole(y_fit,   *popt_dip),  label=f'Dipole fit ($B/y^2$, B={popt_dip[0]:.2f})', color='green')
ax.set_xlabel('Electrode-to-Fiber Distance, $y$ (cm)')
ax.set_ylabel('Peak-to-Peak Amplitude (µV)')
ax.set_title('Peak-to-Peak $\\phi_o$ Amplitude vs Electrode-to-Fiber Distance')
ax.legend()
plt.tight_layout()
plt.savefig('q3_pp_amp_vs_y_cf.png', dpi=300, bbox_inches='tight')
plt.close()

# calculate r^2 values for each model fit
r2_mono = r2_score(pp_amps, monopole(y_values, *popt_mono))
r2_dip  = r2_score(pp_amps, dipole(y_values,   *popt_dip))
print(f'Monopole R²: {r2_mono:.4f}')
print(f'Dipole   R²: {r2_dip:.4f}')