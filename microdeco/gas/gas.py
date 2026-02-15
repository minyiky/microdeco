def _name_from_comp(p_o2: int, p_n2: int, p_he: int) -> str:
    if p_o2 == 21 and p_n2 == 79 and p_he == 0:
        return "Air"
    if p_o2 == 100:
        return "Oxygen"
    if p_he == 0:
        return f"EAN{p_o2}"
    if p_n2 == 0:
        return f"Heliox{p_o2}"
    return f"Tri{p_o2}/{p_he}"


class Gas:
    __slots__ = ('p_o2', 'p_n2', 'p_he', 'name')

    def __init__(self, p_o2: int = 21, p_n2: int = 79, p_he: int = 0):
        if p_o2 + p_n2 + p_he != 100:
            raise ValueError("The sum of the gas percentages must be 100%")
        if p_o2 < 0 or p_n2 < 0 or p_he < 0:
            raise ValueError("Gas percentages must be non-negative")

        self.p_o2 = p_o2 / 100.0
        self.p_n2 = p_n2 / 100.0
        self.p_he = p_he / 100.0

        self.name = _name_from_comp(p_o2, p_n2, p_he)

    def __str__(self) -> str:
        return f"{self.name}"

    def __lt__(self, other):
        if not isinstance(other, Gas):
            return NotImplemented
        # Sort by most oxygen, then least nitrogen
        if self.p_o2 != other.p_o2:
            return self.p_o2 > other.p_o2
        return self.p_n2 < other.p_n2

    def __repr__(self) -> str:
        return f"{self.name} - Gas(Oxygen: {self.p_o2*100:.0f}%, Nitrogen: {self.p_n2*100:.0f}%, Helium: {self.p_he*100:.0f}%)"
