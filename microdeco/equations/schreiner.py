import math

def schreiner(p_alv: float, p_t: float, k: float, t: float, R: float) -> float:
    return p_alv + R * (t - 1 / k) - (p_alv - p_t - (R/k)) * math.exp(-(k*t))
