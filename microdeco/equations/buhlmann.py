def buhlmann(gf: float, p_n2: float, p_he: float, A_n2: float, A_he: float, B_n2: float, B_he: float) -> float:
    p = p_n2 + p_he
    A = (A_n2 * p_n2 + A_he * p_he) / p
    B = (B_n2 * p_n2 + B_he * p_he) / p
    return (p - A * gf) / ((gf / B) + 1 - gf)
