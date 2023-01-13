from .substance import Substance

StoichometricCouple = tuple[float, Substance]

class Reaction:
    reactant: list
    product: list

    def __init__(self, reactant: list[StoichometricCouple], product: list[StoichometricCouple]):
        self.reactant = reactant
        self.product = product
        

    