from diffusion import sim_diffusion
from dolfinx.io.gmshio import model_to_mesh
import gmsh
from mpi4py import MPI
from meshes import generate_2d_rect_mesh, generate_2d_rect_width_holes_mesh

if __name__ == "__main__":
    def generate_full_mesh():
        # https://jsdokken.com/src/tutorial_gmsh.html
        generate_2d_rect_mesh()
        mesh, cell_tags, facet_tags = model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, 2)
        return mesh, cell_tags, facet_tags


    # sim_diffusion(10.0, 1.0, 1.0, generate_full_mesh, 2, 100, "diffusion.xdmf")

    def generate_porous_mesh():
        # https://jsdokken.com/src/tutorial_gmsh.html
        generate_2d_rect_width_holes_mesh(seed= 1000)
        mesh, cell_tags, facet_tags = model_to_mesh(gmsh.model, MPI.COMM_WORLD, 0, 2)
        return mesh, cell_tags, facet_tags

    sim_diffusion(10.0, 1.0, 1.0, generate_porous_mesh, 2, 100, "diffusion_holes.xdmf")