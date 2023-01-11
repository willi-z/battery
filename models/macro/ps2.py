# "pseudo 2d"-model or "Doyle-Fuller-Newmann"-model
from dolfinx import fem, io, mesh, plot
from dolfinx.io.gmshio import model_to_mesh
import gmsh
from mpi4py import MPI
from helpers import Electrode, Electrolyte, Current, generate_ps2_battery_mesh
from mpi4py import MPI
import matplotlib.pyplot as plt


def sim_ps2(
    neg: Electrode,
    sep: Electrolyte,
    pos: Electrode,
    current: Current,
    area: float
    ):
 
    """

    generate_ps2_battery_mesh((16,16))
    mesh, cell_tags, facet_tags = model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, 2)
    fdim = mesh.topology.dim - 1

    s_conc = FiniteElement("CG", mesh.ufl_cell(), 2)
    s_phi = FiniteElement("CG", mesh.ufl_cell(), 2)

    """

    return 0



if __name__ == "__main__":
    omega_neg = mesh.create_rectangle(
            comm=MPI.COMM_WORLD, 
            points=((-1.5, 0.0), (0.5, 1.0)), n=(32, 16),
            cell_type=mesh.CellType.triangle
            )
    omega_sep = mesh.create_interval(
            comm=MPI.COMM_WORLD, 
            points=(-0.5, 0.5), nx=32,
            )
    omega_pos = mesh.create_rectangle(
            comm=MPI.COMM_WORLD, 
            points=((0.5, 0.0), (1.5, 1.0)), n=(32, 16),
            cell_type=mesh.CellType.triangle
            )
    
    omega_all = omega_neg + omega_sep + omega_pos

    plt.plot(omega_all)
    plt.show()
    plt.savefig('omega.png')