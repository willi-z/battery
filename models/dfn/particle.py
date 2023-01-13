import numpy as np
import multiprocessing as mp


class Particle:
    r_particle: float       # particle radius [m]
    diffusivity: float      # solid diffusivity [m^2/s]
    areas_shell : np.ndarray
    volumes_shell : np.ndarray
    dr: float

    def __init__(self, r_particle: float, diffusivity: float, n_shells: float):
        self.r_particle = r_particle
        self.diffusivity = diffusivity
        self.dr = r_particle / n_shells
        self.areas_shell = 4 * np.pi * (self.dr * np.arange(1, n_shells + 1, 1)) ** 2
        self.volumes_shell = (4/3.) * np.pi * ((self.dr * np.arange(1, n_shells + 1, 1)) ** 3 - (self.dr * np.arange(0, n_shells, 1)) ** 3)

    def get_n_shells(self) -> int:
        return len(self.areas_shell)

    def get_radius(self) -> int:
        return self.r_particle


class ParticleSolver:
    particle: Particle
    # c: np.ndarray # concentration c[0] concentration at center, c[-1] concntration at surface

    def __init__(self, particle: Particle):
        self.particle = particle

    def simulate_dt(self, c: np.ndarray, j_ext: float, dt: float) -> float:
        """
        Returns surface concentration
        """
        diffusivity = self.particle.diffusivity
        areas_shell = self.particle.areas_shell
        volumes_shell = self.particle.volumes_shell

        N = -diffusivity * np.diff(c) / self.particle.dr  # flux density at surface between "shells"
        M = N * areas_shell[:-1]           # total moles crossing surface between "shells"
        c = c + (np.append([0], M) - np.append(M, [0])) * dt/volumes_shell
        c[-1] = c[-1] - j_ext * areas_shell[-1] * dt/volumes_shell[-1]   # boundary
        return c[-1]



class ParallelParticleSolver:
    solver: ParticleSolver
    cs = np.array
    time: np.ndarray

    def __init__(self, particle: Particle, c_0: float, dofs: list):
        self.solver = ParticleSolver(particle)
        n_solver = dofs # len(dofs)
        self.cs = np.ones((n_solver, particle.get_n_shells())) * c_0
        self.time = np.ones(n_solver)
        

    def simulate_dt(self, j_ext: np.ndarray, dt: float) -> np.ndarray:
        inputs =  map(lambda x0, x1, x2: (x0, x1, x2), self.cs, j_ext, self.time * dt)
        with mp.Pool(processes=4) as pool:
            results = pool.starmap(self.solver.simulate_dt, inputs)
        return results

if __name__ == "__main__":
    particle = Particle(10e-6, 1e-14, 10)
    c_0 = 9600          # initial Li concentration [mol/m^3]
    solver = ParallelParticleSolver(particle, c_0, 5)
    j_0 = 5000*particle.get_radius()/3/1800 # Li flux [mol/m^2/s]
    j_ext = (np.ones(10) + np.ndarray(10)/ 10.0) * j_0
    print(solver.cs)
    result = solver.simulate_dt(j_ext, 1)
    print(solver.cs)
    print(result)
