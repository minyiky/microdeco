class GradientFactors:

    def __init__(self, gf_low: float = 0.3, gf_high: float = 0.85):
        if not (0 < gf_low <= gf_high <= 1):
            raise ValueError("Gradient factors must satisfy 0 < gf_low <= gf_high <= 1")
        self.gf_low = gf_low
        self.gf_high = gf_high
        self.initial_depth = 0.0

    def set_initial_depth(self, initial_depth: float) -> None:
        self.initial_depth = initial_depth

    def get_current_gf(self, current_depth: float) -> float:
        if current_depth > self.initial_depth:
            return self.gf_low
        ratio = (self.initial_depth - current_depth) / self.initial_depth
        return self.gf_low + (self.gf_high - self.gf_low) * ratio
