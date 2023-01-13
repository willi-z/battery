from .electrode import Electrode
from .electrolyte import Electrolyte
from .reaction import Reaction


class Battery:
    neg: Electrode
    pos: Electrode
    sep: Electrolyte
    react: Reaction

    def __init__(self, neg, pos, sep, react):
        pass

    