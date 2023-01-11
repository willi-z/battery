import numpy as np

from dolfinx import fem, io
from ufl import ds, dx, grad, dot, inner, TrialFunction, TestFunction, FiniteElement, Measure, FacetNormal, lhs, rhs, SpatialCoordinate
from petsc4py import PETSc
from petsc4py.PETSc import ScalarType
from tqdm import tqdm


def get_vals_at_tag() -> np.ndarray:
    pass

def sim_diffusion(
    diff_const: float,
    eps_const: float,
    conc_inflow_const: float,
    mesh_generator,
    model_rank,
    T_end,
    fname: str,
    num_steps = 100,
    substance_tag = 1,
    inflow_tag = 2,
    outflow_tag = 3
    ):
    """
    simulate eps*dc/dt+ \nabla \cdot (-diffusivity cdot \nabla c) = 
    """
    
    mesh, cell_tags, facet_tags = mesh_generator()

    fdim = mesh.topology.dim - 1
    s_conc = FiniteElement("CG", mesh.ufl_cell(), 2) # scaler valued concentration function
    fs_conc = fem.FunctionSpace(mesh, s_conc) # functions space concentration

    x = SpatialCoordinate(mesh)

    eps = fem.Constant(mesh, ScalarType(eps_const))
    diff = fem.Constant(mesh, ScalarType(diff_const))

    conc_inflow = fem.Function(fs_conc)
    conc_inflow_func = lambda x: 1.0 * np.ones(len(x[0]))
    conc_inflow.interpolate(conc_inflow_func)
    bc_conc_inflow = fem.dirichletbc(conc_inflow, fem.locate_dofs_topological(fs_conc, fdim, facet_tags.find(inflow_tag)))
    bcs = [bc_conc_inflow]

    # Solver
    # ds = Measure("ds", domain=e, subdomain_data=facet_tag)
    t = 0
    T = T_end
    num_steps = 100
    dt = (T - t) / num_steps

    conc_n = fem.Function(fs_conc)
    conc_n.name = "conc"

    conc_old = fem.Function(fs_conc)
    conc_old.name = "conc_old"

    conc, v = TrialFunction(fs_conc), TestFunction(fs_conc) # trail and test

    a = eps* conc * v *dx + dt * dot(diff * grad(conc), grad(v)) * dx
    L = eps * inner(conc_n, v) * dx

    bilinear_form = fem.form(a)
    linear_form = fem.form(L)

    A = fem.petsc.assemble_matrix(bilinear_form, bcs=bcs)
    A.assemble()
    b = fem.petsc.create_vector(linear_form)


    solver = PETSc.KSP().create(mesh.comm)
    solver.setOperators(A)
    solver.setType(PETSc.KSP.Type.PREONLY)
    solver.getPC().setType(PETSc.PC.Type.LU)

    xdmf = io.XDMFFile(mesh.comm, fname, "w")
    xdmf.write_mesh(mesh)
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
