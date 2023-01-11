import gmsh
import numpy as np

def generate_ps2_battery_mesh(n_difs: tuple):
    gmsh.initialize()
    gmsh.clear()

    gmsh.model.add("ps2 battery mesh") # title

    background = gmsh.model.occ.addRectangle(0, 0, 0, 3.0, 1.0)

    x, y, z = 0, 0, 0
    neg = gmsh.model.occ.addRectangle(x, y, z, 1.0, 1.0, 2)

    x = 1.0
    sep = gmsh.model.occ.addRectangle(x, y, z, 1.0, 1.0, 3)

    x = 2.0
    pos = gmsh.model.occ.addRectangle(x, y, z, 1.0, 1.0, 4)

    gdim = 2
    all_surfaces = [(gdim, neg), (gdim, sep), (gdim, pos)]
    whole_domain = gmsh.model.occ.fragment([(2, background)], all_surfaces)
    # battery = gmsh.model.occ.merge([(gdim, sep)], [(gdim, neg), (gdim, pos)])
    gmsh.model.occ.synchronize()

    volumes = gmsh.model.getEntities(dim=gdim)
    print(len(volumes))
    tag_neg = 1
    gmsh.model.addPhysicalGroup(volumes[0][0], [volumes[0][1]], tag= tag_neg)
    gmsh.model.setPhysicalName(volumes[0][0], tag_neg, "neg")

    tag_sep = 2
    gmsh.model.addPhysicalGroup(volumes[1][0], [volumes[1][1]], tag= tag_sep)
    gmsh.model.setPhysicalName(volumes[1][0], tag_sep, "sep")

    tag_pos = 3
    gmsh.model.addPhysicalGroup(volumes[2][0], [volumes[2][1]], tag= tag_pos)
    gmsh.model.setPhysicalName(volumes[2][0], tag_pos, "pos")
    # assert(len(volumes) == 1)

    bound_neg, int_neg_sep, int_sep_pos, bound_pos = [], [], [], []
    tag_bound_neg, tag_int_neg_sep, tag_int_sep_pos, tag_bound_pos = 4, 5, 6, 7


    surfaces = gmsh.model.getEntities(dim = 2)
    boundaries = gmsh.model.getBoundary(surfaces, oriented=False)
    for boundary in boundaries:
        center_of_mass = gmsh.model.occ.getCenterOfMass(boundary[0], boundary[1])
        if np.allclose(center_of_mass[0], 0):
            bound_neg.append(boundary[1])
        elif np.allclose(center_of_mass[0], 1.0):
            int_neg_sep.append(boundary[1])
        elif np.allclose(center_of_mass[0], 2.0):
            int_sep_pos.append(boundary[1])
        elif np.allclose(center_of_mass[0], 3.0):
            bound_pos.append(boundary[1])

    gmsh.model.addPhysicalGroup(1, bound_neg, tag=tag_bound_neg, name="boundary neg")
    gmsh.model.addPhysicalGroup(1, int_neg_sep, tag=tag_int_neg_sep, name="interface neg-sep")
    gmsh.model.addPhysicalGroup(1, int_sep_pos, tag=tag_int_sep_pos, name="interface sep-pos")
    gmsh.model.addPhysicalGroup(1, bound_pos, tag=tag_bound_pos, name="boundary pos")
    

    gmsh.model.mesh.generate(2) # 2D mesh
    gmsh.model.mesh.setOrder(1) # linear elments
    gmsh.model.mesh.refine()


if __name__=="__main__":
    generate_ps2_battery_mesh((20,10))
    gmsh.write("ps2.msh")