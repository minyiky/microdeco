'''
'''

from microdeco.equations import buhlmann

class Tissue:
    def __init__(self, H_n2: float, A_n2: float, B_n2: float) -> None:
        self.H_n2, self.A_n2, self.B_n2 = H_n2, A_n2, B_n2

    def buhlmann(self, gf: float, pp_n2: float, pp_he: float) -> float:
        return buhlmann(gf, pp_n2, pp_he, self.A_n2, self.B_n2, 0, 0)


class Buhlmann_ZH16_C:
    H_n2 = [5.0, 8.0, 12.5, 18.5, 27.0, 38.3, 54.3, 77.0, 109.0,
            146.0, 187.0, 239.0, 305.0, 390.0, 498.0, 635.0]
    A_n2 = [1.1696, 1.0000, 0.8618, 0.7562, 0.6200, 0.5043, 0.4410, 0.4000,
            0.3750, 0.3500, 0.3295, 0.3065, 0.2835, 0.2610, 0.2480, 0.2327]
    B_n2 = [0.5578, 0.6514, 0.7222, 0.7825, 0.8126, 0.8434, 0.8693, 0.8910,
            0.9092, 0.9222, 0.9319, 0.9403, 0.9477, 0.9544, 0.9602, 0.9653]

    tissues = [
        Tissue(5.0, 1.1696, 0.5578),
        Tissue(8.0, 1.0000, 0.6514),
        Tissue(12.5, 0.8618, 0.7222),
        Tissue(18.5, 0.7562, 0.7825),
        Tissue(27.0, 0.6200, 0.8126),
        Tissue(38.3, 0.5043, 0.8434),
        Tissue(54.3, 0.4410, 0.8693),
        Tissue(77.0, 0.4000, 0.8910),
        Tissue(109.0, 0.3750, 0.9092),
        Tissue(146.0, 0.3500, 0.9222),
        Tissue(187.0, 0.3295, 0.9319),
        Tissue(239.0, 0.3065, 0.9403),
        Tissue(305.0, 0.2835, 0.9477),
        Tissue(390.0, 0.2610, 0.9544),
        Tissue(498.0, 0.2480, 0.9602),
        Tissue(635.0, 0.2327, 0.9653),
    ]
    
    def ceil_list(self) -> float:
         return max(
            buhlmann(1, 0.78, 0, self.A_n2[i], self.B_n2[i], 0, 0)
            for i in range(16)
        )

    def ceil_class(self) -> float:
         return max(
            buhlmann(1, 0.78, 0, self.tissues[i].A_n2, self.tissues[i].B_n2, 0, 0)
            for i in range(16)
        )
    
    def ceil_class_func(self) -> float:
         return max(
            self.tissues[i].buhlmann(1, 0.78, 0)
            for i in range(16)
        )
    
