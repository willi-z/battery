from skfem import *
from skfem.helpers import dot, grad
import numpy as np
from tqdm import tqdm


from os.path import splitext
from sys import argv
from skfem.visuals.matplotlib import draw, plot


# PARAMETER
r_particle = 10e-6  # particle radius [m]

c_max = 12000       # maximum Li concentration [mol/m^3]
c_0 = 9600          # initial Li concentration [mol/m^3]
assert c_0 <= c_max

j_0 = 5000*r_particle/3/1800 # Li flux [mol/m^2/s]
diff_particle = 1e-14       # solid diffusivity [m^2/s]

dt = 1                      # time steps [s]
t_discharge = 0.5 * 60 * 60 # simulation time [s]
t_rest = 1.0 * 60 * 60      # simulation time [s]
t_charge = 0.5 * 60 * 60    # simulation time [s]

js_discharge = np.append(j_0 * np.ones(int(t_discharge)) , np.zeros(int(t_rest)))
js_charge = np.append(-j_0 * np.ones(int(t_charge)) , np.zeros(int(t_rest)))

js_external = np.append(js_discharge, js_charge)
n_steps = len(js_external)


# SIMULATION CONTROL
n_points_particle = 10 +1

# MESHING
space = np.linspace(0, 1, n_points_particle)
mesh = (MeshLine(space)
      .with_boundaries({
         'inner': lambda xi: xi[0] == 0,
         'surface': lambda xi: xi[0] == 1.0
      })
   )

basis = Basis(mesh, ElementLineP1())
basis_surf = basis.boundary('surface')

# ASSEMBLE A
@BilinearForm
def diffusion(u, v, w):
   r = w.x[0]
   return r**2 * r_particle * u * v + \
            dt * diff_particle * r**2/r_particle * dot(grad(u), grad(v))

A = asm(diffusion, basis)

@BilinearForm
def L(conc, v, w):
     r = w.x[0]
     return r**2 * r_particle * v * conc

@LinearForm
def L_surf(v, w):
     return dt * dot(w.n, v) * r_particle #  * -js_external[i]

M = L.assemble(basis)
f_surf = L_surf.assemble(basis_surf)

dofs = basis.get_dofs("surface")
# basis_surf = FacetBasis(mesh, basis.elem).


# INTIAL CONDITIONS
def inital(x):
   return np.ones((x.shape[1], x.shape[2])) * c_0

c = basis.project(inital)

c_surface = np.zeros(n_steps) # array of concentration at surface for each timestep
c_surface[0] = c_0


for i in tqdm(range(n_steps - 1)):
   f_1 = M.dot(c)
   f_2 = f_surf * -js_external[i]
   c_n = solve(A, f_1 + f_2)
   c = c_n
   c_surf = c[basis.get_dofs("surface")]
   c_surface[i+1] = c_surf[0]


import matplotlib.pyplot as plt
t = np.arange(n_steps)
plt.plot(t/(60*60), c_surface)
# plt.show()
plt.savefig('spm_fast.png')
