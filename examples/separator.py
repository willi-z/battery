from dolfinx import fem, io, mesh, plot
from mpi4py import MPI
import numpy as np
from ufl import ds, dx, grad, dot, inner, TrialFunction, TestFunction, Measure, FacetNormal, lhs, rhs, SpatialCoordinate
from petsc4py import PETSc
from petsc4py.PETSc import ScalarType
from tqdm import tqdm


x_a = 0.0
x_b = 2.0

t = 0
T = 10
num_steps = 100
dt = (T - t) / num_steps

msh = mesh.create_rectangle(comm=MPI.COMM_WORLD,
                            points=((x_a, 0.0), (x_b, 1.0)), n=(32, 16),
                            cell_type=mesh.CellType.triangle,)
V = fem.FunctionSpace(msh, ("Lagrange", 1))
x = SpatialCoordinate(msh)
conc, v = TrialFunction(V), TestFunction(V)


eps = fem.Constant(msh, ScalarType(0.5))
diffusion = fem.Constant(msh, ScalarType(1))

def initial_condition(x):
    return 0.0 * x[0]

conc_n = fem.Function(V)
conc_n.name = "conc"
conc_n.interpolate(initial_condition)

conc_old = fem.Function(V)
conc_old.name = "conc_old"
conc_old.interpolate(initial_condition)

# BC's
boundaries = [(1, lambda x: np.isclose(x[0], x_a)),
              (2, lambda x: np.isclose(x[0], x_b))]


facet_indices, facet_markers = [], []
fdim = msh.topology.dim - 1
for (marker, locator) in boundaries:
    facets = mesh.locate_entities(msh, fdim, locator)
    facet_indices.append(facets)
    facet_markers.append(np.full_like(facets, marker))
facet_indices = np.hstack(facet_indices).astype(np.int32)
facet_markers = np.hstack(facet_markers).astype(np.int32)
sorted_facets = np.argsort(facet_indices)
facet_tag = mesh.meshtags(msh, fdim, facet_indices[sorted_facets], facet_markers[sorted_facets])


ds = Measure("ds", domain=msh, subdomain_data=facet_tag)

a = eps* conc * v *dx + dt * dot(grad(conc), grad(v)) * dx
L = eps * inner(conc_n, v) * dx

class BoundaryCondition():
    def __init__(self, type, marker, values):
        self._type = type
        if type == "Dirichlet":
            u_D = fem.Function(V)
            u_D.interpolate(values)
            facets = facet_tag.find(marker)
            dofs = fem.locate_dofs_topological(V, fdim, facets)
            self._bc = fem.dirichletbc(u_D, dofs)
        elif type == "Neumann":
                self._bc = inner(values, v) * ds(marker)
        else:
            raise TypeError("Unknown boundary condition: {0:s}".format(type))
    @property
    def bc(self):
        return self._bc

    @property
    def type(self):
        return self._type

conc_dc = lambda x: 1.0 * np.ones(len(x[0]))
val_nc = lambda x: 0.0 * x[0]

boundary_conditions = [
    BoundaryCondition("Dirichlet", 1, conc_dc),
    ]

bcs = []
for condition in boundary_conditions:
    if condition.type == "Dirichlet":
        bcs.append(condition.bc)
    else:
        L += condition.bc


bilinear_form = fem.form(a)
linear_form = fem.form(L)

A = fem.petsc.assemble_matrix(bilinear_form, bcs=bcs)
A.assemble()
b = fem.petsc.create_vector(linear_form)


solver = PETSc.KSP().create(msh.comm)
solver.setOperators(A)
solver.setType(PETSc.KSP.Type.PREONLY)
solver.getPC().setType(PETSc.PC.Type.LU)

xdmf = io.XDMFFile(msh.comm, "diffusion.xdmf", "w")
xdmf.write_mesh(msh)
xdmf.write_function(conc_n, t)

print(bcs)

for i in tqdm(range(num_steps)):
    t += dt

    # Update the right hand side reusing the initial vector
    with b.localForm() as loc_b:
        loc_b.set(0)
    fem.petsc.assemble_vector(b, linear_form)
    
    # Apply Dirichlet boundary condition to the vector
    fem.petsc.apply_lifting(b, [bilinear_form], bcs = [bcs])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
    fem.petsc.set_bc(b, bcs)

    # Solve linear problem
    solver.solve(b, conc_old.vector)
    conc_old.x.scatter_forward()

    # Update solution at previous time step (u_n)
    conc_n.x.array[:] = conc_old.x.array

    # Write solution to file
    xdmf.write_function(conc_n, t)


xdmf.close()