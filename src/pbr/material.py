"""(5) Material balance per tube: dN_x/dz = -A_c r   <=>   dX/dz = A_c r / N_x0.
(r >= 0 is the o-xylene consumption rate per unit bed volume; X = 1 - N_x/N_x0.)"""


def coeff(p):
    """A_c / N_x0 [m2 s/mol]."""
    return p.A_c / p.N_x0


def dXdz(r, p):
    return coeff(p) * r
