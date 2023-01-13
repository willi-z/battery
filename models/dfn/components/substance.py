

class Substance:
    name: str  # [Optional]
    charge: int  # +: 1, -: -1
    def __init__(charge: int, name = ''):
        self.name = name
        self.charge = charge
