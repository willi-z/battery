import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm


"""
Simulate the diffusion of lithium within a single solid particle.
Here a finite volume method for discretizising the diffusion equation is used.
A single solid particle is divided into spherical shells having equal thickness (like an onion).
At each time step the total flux of lithium from  one shell to another is calculated,
and the concentration of lithium within each sell is then updated.
Lithium is forced into or out of the outermost shell via an imposed cell current.

The partical is discharged, rested, charged, and rested over 3 hours.

Original code from "Battery Modeling" (Volume 1, p. 128) by Gregory L. Plett
"""

r_particle = 10e-6  # particle radius [m]

c_max = 12000       # maximum Li concentration [mol/m^3]
c_0 = 9600          # initial Li concentration [mol/m^3]
assert c_0 <= c_max

j_0 = 5000*r_particle/3/1800 # Li flux [mol/m^2/s]
diffusivity = 1e-14   # solid diffusivity [m^2/s]

dt = 1                      # time steps [s]
t_discharge = 0.5 * 60 * 60 # simulation time [s]
t_rest = 1.0 * 60 * 60      # simulation time [s]
t_charge = 0.5 * 60 * 60    # simulation time [s]

js_discharge = np.append(j_0 * np.ones(int(t_discharge)) , np.zeros(int(t_rest)))
js_charge = np.append(-j_0 * np.ones(int(t_charge)) , np.zeros(int(t_rest)))

js_external = np.append(js_discharge, js_charge)
n_steps = len(js_external)

# Simulation control

n_shells = 20   # number of "shells radially"
dr = r_particle / n_shells # width of each "shell"
areas_surface = 4 * np.pi * (dr * np.arange(1, n_shells + 1, 1)) ** 2
volumes_shell = (4/3.) * np.pi * ((dr * np.arange(1, n_shells + 1, 1)) ** 3 - (dr * np.arange(0, n_shells, 1)) ** 3)


c = c_0 * np.ones(n_shells) # concentration profile versius "r" domain
c_surface = np.zeros(n_steps) # array of concentration at surface for each timestep
c_surface[0] = c_0

for i in tqdm(range(n_steps - 1)):
    N = -diffusivity * np.diff(c) / dr  # flux density at surface between "shells"
    M = N * areas_surface[:-1]           # total moles crossing surface between "shells"
    c = c + (np.append([0], M) - np.append(M, [0])) * dt/volumes_shell
    c[-1] = c[-1] - js_external[i] * areas_surface[-1] * dt/volumes_shell[-1]   # boundary
    c_surface[i+1] = c[-1]


x = c_surface / c_max

u_ocp = lambda x: 3.7 * np.ones(len(x))

t = np.arange(n_steps)
plt.plot(t/(60*60), c_surface)
# plt.show()
plt.savefig('spm.png')