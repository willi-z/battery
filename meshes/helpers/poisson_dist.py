import numpy as np
import random

import matplotlib.pyplot as plt
# https://sighack.com/post/poisson-disk-sampling-bridsons-algorithm


class Grid():
    array: list

    # static
    sizes: list
    offsets: list 
    length: int
    cell_size: float
    dim: int
    combis: np.ndarray
    periodic: bool
    def __init__(self, sizes: list, cell_size: float, perodic=False) -> None:
        self.dim = len(sizes)
        self.sizes = np.array(sizes, dtype=np.uintc)
        self.cell_size = cell_size
        self.periodic = perodic
        self.length = 1
        self.offsets = [1] * self.dim
        # print(sizes)
        for i in range(self.dim):
            self.length *= sizes[i]
            for j in range(i + 1, self.dim):
                self.offsets[j] *= self.sizes[i]
                # print(f"{j} | {i} = {self.offsets[j]}")
        # print(self.offsets)
        self.array = [None] * int(self.length)
        self.combis = get_combinations([-2, -1, 0, 1, 2], self.dim)

        corners = []
        for i in range(len(self.combis)):
            if np.isclose(np.sum(np.square(self.combis[i])), np.square(2) * self.dim):
                corners.append(i)
        self.combis = np.delete(self.combis, corners, 0)

    def get_id(self, indeces: list) -> int:
        assert len(indeces) == self.dim
        for i in range(self.dim):
            if indeces[i] >= self.sizes[i] or indeces[i] < 0:
                return None
        index = 0
        for i in range(self.dim):
            index += int(indeces[i]) * self.offsets[i] 
        return index

    def calc_indeces(self, point: np.ndarray) -> list:
        assert len(point) == self.dim
        return np.floor(point / self.cell_size)

    def get_cell(self, indeces: list) -> np.ndarray:
        id = self.get_id(indeces)
        if id is None:
            return None
        return self.array[id]


    def insert_point(self, point: np.ndarray):
        id = self.get_id(self.calc_indeces(point))
        assert id is not None
        self.array[id] = point

    def get_neighbors(self, point: np.ndarray):
        assert len(point) == self.dim
        indeces = np.floor(point / self.cell_size)
        neighbors = []
        for combi in self.combis:
            indeces_n = indeces + combi
            if self.periodic:
                for i in range(len(indeces_n)):
                    if indeces[i] < 0:
                        indeces[i] += sizes[i]
                    if indeces[i] >= sizes[i]:
                        indeces[i] -= sizes[i]
            neighbor = self.get_cell(indeces_n)
            if neighbor is not None:
                neighbors.append(neighbor)
        return neighbors

    def is_valid(self, point: np.ndarray, min_distance: float) -> bool:
        assert len(point) == self.dim
        for i in range(self.dim):
            if point[i] < 0 or point[i] > self.sizes[i] * self.cell_size:
                return False
        neighbors = self.get_neighbors(point)
        radius_2 = min_distance ** 2
        for neighbor in neighbors:
            if  np.sum(np.square(neighbor - point)) < radius_2:
                return False
        # for neighbor in neighbors:
        #     radius = np.sum(np.square(neighbor - point))
        #     print(f"{radius_2} | {radius} : {neighbor} | {point}")
        return True


def get_combinations(vals: list, dim: int) -> list:
    length = len(vals) ** dim
    combinations = np.zeros((length, dim))
    iterator = np.zeros(dim, dtype=np.uintc)
    for i in range(length):
        for j in range(dim):
            combinations[i][j] = vals[iterator[j]]
        c = 0
        while True:
            iterator[c] += 1
            if iterator[c] < len(vals):
                break
            else:
                iterator[c] = 0
                c += 1
                if c >= dim:
                    break
    return combinations

def poision_distribution(limits: list, min_dist: float, max_tries = 10) -> list:
    dim = len(limits)
    cell_size = min_dist / np.sqrt(dim) # cell size so each grid contains at most one sample
    sizes = np.ceil(np.array(limits) / cell_size) + 1 # array with the number of cells per dimension

    grid = Grid(sizes, cell_size)
    point_0 = np.random.rand(dim) * np.array(limits)

    active = []
    points = []
    grid.insert_point(point_0)
    active.append(point_0)
    points.append(point_0)

    while len(active) > 0:
        n_active = len(active)
        choice = random.randrange(n_active)
        point = active[choice]
        
        found = False
        for i in range(max_tries):
            theta = random.random() * np.pi
            phi = random.random() * 2 * np.pi
            radius = random.random() * min_dist + min_dist
            if dim == 1:
                theta = 0.0
            p_new = [radius * np.cos(theta)]
            if dim > 1: 
                p_new.append(radius * np.sin(theta) * np.cos(phi))
            if dim == 3:
                p_new.append(radius * np.sin(theta) * np.sin(phi))
            point_new = np.array(p_new) + point

            if grid.is_valid(point_new, min_dist):
                points.append(point_new)
                # print(len(points))
                grid.insert_point(point_new)
                active.append(point_new)
                found = True
                break
        if not found:
            del active[choice]
    return points


def test_detection():
    min_dist = 5
    limits = [10, 10]
    dim = len(limits)
    cell_size = min_dist / np.sqrt(dim)
    sizes = np.ceil(np.array(limits) / cell_size) + 1 # array with the number of cells per dimension

    grid = Grid(sizes, cell_size)

    points = [
        np.array([4.69984105, 3.62739606]), 
        np.array([9.96902699, 3.91896774]), # 1
        np.array([6.02738207, 8.46495918]), # 2
        np.array([10.76294467,  9.43061412]), # 3
        np.array([ 0.38163462, 13.08034709]), # 4
        np.array([ 5.30152173, 13.99382399]) # 5
    ]
    print(f"teste {len(points)}")
    print(np.sum(np.square(points[3] - points[2])) - min_dist**2)
    grid.insert_point(points[0])
    for i in range(1, len(points)):
        if not grid.is_valid(points[i], min_dist):
            print(f"{i}: {points[i]}")
        grid.insert_point(points[i])
    return points

if __name__ == "__main__":
    # print(get_combinations([-1, 0, 1], 2))
    radius_min = 1
    limits = [10, 10]
    points = poision_distribution(limits, radius_min, 30)
    print(points)
    # points = test_detection()
    x = []
    y = []
    for point in points:
        x.append(point[0])
        y.append(point[1])
        circle = plt.Circle((point[0], point[1]), radius_min/2, fill = False)
        plt.gcf().gca().add_artist(circle)
    plt.scatter(x, y)

    cell_size = radius_min / np.sqrt(2)
    sizes = np.ceil(np.array(limits) / cell_size) + 1
    plt.xticks(np.arange(sizes[0])*cell_size)
    plt.yticks(np.arange(sizes[1])*cell_size)
    plt.grid(True)

    plt.show()
    plt.savefig('possion_dist.png')
