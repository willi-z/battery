# https://scikit-fem.readthedocs.io/en/latest/
import skfem as fem
from skfem.helpers import dot, grad  # helpers make forms look nice

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

n_r = 20 + 1   # number of "shells radially"
dr = r_particle / n_r # width of each "shell"

space = np.linspace(0, 1, n_r)
mesh = (fem.MeshLine(space)
    .with_boundaries({
        'center': lambda xi: xi[0] == 0.,
        'surface': lambda xi: xi[0] == 1.0,
    })
    )
basis = fem.Basis(mesh, fem.ElementLineP1())

@fem.BilinearForm
def a(u, v, w):
    # u = c
    r = w.x[0] # global coordinates
    return r**2 * r_particle**2 * u * v + \
        dt* diffusivity * r**2/r_particle * dot(grad(u), grad(v))

A = a.assemble(basis)

def inital(x):
    return np.ones(x.shape) * c_0

c = basis.project(inital)

# c = np.concatenate((C0, np.zeros(basis.N)))

c_surface = np.zeros(n_steps) # array of concentration at surface for each timestep
c_surface[0] = c_0

for i in tqdm(range(n_steps - 1)):
    conc = basis.interpolate(c)
    @fem.LinearForm
    def L(v, w):
        r = w.x[0]  # global coordinates
        f = r**2 * r_particle**2 * conc
        return f * v

    @fem.LinearForm
    def L_neumann(v, w):
        return -js_external[i] * dot(w.n, v)

    l = L.assemble(basis) + L_neumann.assemble(basis[])
    # cond = fem.condense(A, l, D=D)
    c = fem.solve(A, l)
    c_surf = c[basis.get_dofs("surface")]
    c_surface[i+1] = c_surf


x = c_surface / c_max

u_ocp = lambda x: 3.7 * np.ones(len(x))

t = np.arange(n_steps)
plt.plot(t/(60*60), c_surface)
plt.show()
plt.savefig('spm_fem.png')