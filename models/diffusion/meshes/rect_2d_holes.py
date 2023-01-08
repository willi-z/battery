import gmsh
import numpy as np

import sys
from helpers import poision_distribution

def generate_2d_plane_width_holes_mesh(
    width = 30, 
    height = 10, 
    r_holes = 1.5, 
    min_dist=0.5,
    seed=None):
    gmsh.initialize()
    gmsh.clear()

    gmsh.model.add("2d plane") # title

    x, y, z = 0, -height/2, 0
    plane = gmsh.model.occ.addRectangle(x, y, z, width, height, 0)
    obstacles = []

    """
    References:
        Supersampling 
        https://en.wikipedia.org/wiki/Supersampling

        Poisson Disk Sampling + Perlin Noise
        http://devmag.org.za/2009/05/03/poisson-disk-sampling/
    """

    points = poision_distribution([width, height], 2 * r_holes + min_dist, 20, seed)
    n_holes = len(points)
    obstacles = []
    cuts = []
    gdim = 2
    offset = - np.array([0, height/2])
    for i in range(n_holes):
        center = points[i] + offset
        obstacle = gmsh.model.occ.addDisk(center[0], center[1], 0, r_holes, r_holes)
        cuts.append((gdim, obstacle))
        obstacles.append(obstacle)
    substance = gmsh.model.occ.cut([(gdim, plane)], cuts)
    gmsh.model.occ.synchronize()

    """
    0 - point
    1 - line
    2 - surface
    3 - volume
    """
    substance_marker = 1
    volumes = gmsh.model.getEntities(dim=gdim)
    assert(len(volumes) == 1)
    gmsh.model.addPhysicalGroup(volumes[0][0], [volumes[0][1]], tag= substance_marker)
    gmsh.model.setPhysicalName(volumes[0][0], substance_marker, "substance")
    # gmsh.model.addPhysicalGroup(2, [plane], tag=1, name="rect plane")

    
    inflow, outflow, walls, obstacles = [], [], [], []
    inflow_tag, outflow_tag, wall_tag, obstacle_tag = 2, 3, 4, 5

    surfaces = gmsh.model.getEntities(dim = 2)
    assert len(surfaces) == 1
    boundaries = gmsh.model.getBoundary(surfaces, oriented=False)
    for boundary in boundaries:
        center_of_mass = gmsh.model.occ.getCenterOfMass(boundary[0], boundary[1])
        if np.allclose(center_of_mass[0], 0):
            inflow.append(boundary[1])
        elif np.allclose(center_of_mass[0], width):
            outflow.append(boundary[1])
        elif np.allclose(center_of_mass[1], -height/2) or np.allclose(center_of_mass[1], height/2):
            walls.append(boundary[1])
        else:
            obstacles.append(boundary[1])

    gmsh.model.addPhysicalGroup(1, inflow, tag=inflow_tag, name="inflow")
    gmsh.model.addPhysicalGroup(1, outflow, tag=outflow_tag, name="outflow")
    gmsh.model.addPhysicalGroup(1, walls, tag=wall_tag, name="walls")
    gmsh.model.addPhysicalGroup(1, obstacles, tag=obstacle_tag, name="obstacle")
    
    gmsh.model.mesh.generate(2) # 2D mesh
    gmsh.model.mesh.setOrder(1) # linear elments
    gmsh.model.mesh.refine()

if __name__=="__main__":
    generate_2d_plane_width_holes_mesh()
    gmsh.write("2d_plane+holes.msh")