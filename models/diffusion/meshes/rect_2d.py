import gmsh
import numpy as np


def generate_2d_plane_mesh(width = 30, height = 10):
    gmsh.initialize()
    gmsh.clear()

    gmsh.model.add("2d plane") # title

    x, y, z = 0, -height/2, 0
    plane = gmsh.model.occ.addRectangle(x, y, z, width, height, 0)
    gmsh.model.occ.synchronize()

    """
    0 - point
    1 - line
    2 - surface
    3 - volume
    """
    gmsh.model.addPhysicalGroup(2, [plane], tag=1, name="substance")

    inflow, outflow, walls = [], [], []
    inflow_tag, outflow_tag, wall_tag = 2, 3, 4

    surfaces = gmsh.model.getEntities(dim = 2)
    assert len(surfaces) == 1
    boundaries = gmsh.model.getBoundary(surfaces, oriented=False)
    for boundary in boundaries:
        center_of_mass = gmsh.model.occ.getCenterOfMass(boundary[0], boundary[1])
        if np.allclose(center_of_mass, [0, 0, 0]):
            inflow.append(boundary[1])
        elif np.allclose(center_of_mass, [width, 0, 0]):
            outflow.append(boundary[1])
        elif np.allclose(center_of_mass, [width/2, -height/2, 0]) or np.allclose(center_of_mass, [width/2, height/2, 0]):
            walls.append(boundary[1])

    gmsh.model.addPhysicalGroup(1, inflow, tag=inflow_tag, name="inflow")
    gmsh.model.addPhysicalGroup(1, outflow, tag=outflow_tag, name="outflow")
    gmsh.model.addPhysicalGroup(1, walls, tag=wall_tag, name="walls")

    gmsh.model.mesh.generate(2) # 2D mesh
    gmsh.model.mesh.setOrder(1) # linear elments
    gmsh.model.mesh.refine()


if __name__=="__main__":
    generate_2d_plane_mesh()
    gmsh.write("2d_plane.msh")