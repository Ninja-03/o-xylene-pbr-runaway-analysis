"""(6) Energy balance per tube (constant c_p, mass basis, salt at constant T_c):

    Qrho c_p dT/dz = (-dH) A_c r - 2 pi R_tube U (T - T_c)

Heat-release definition: q_rx = (-dH) A_c r [W/m] (> 0, exothermic); total = (-dH) N_x0 X(L).
Cooling: q_cool = 2 pi R_tube U (T - T_c) [W/m]. SIGN CONVENTION: q_cool > 0 means heat flows from the
gas to the salt (removed from the reacting gas); Q_cool = int_0^L q_cool dz is the cooling duty.
Global form: Qrho c_p [T(L) - T0] = (-dH) N_x0 X(L) - Q_cool.
"""
from . import properties


def heat_capacity_flow(p):
    """Qrho c_p [W/K]."""
    return p.Qrho * properties.cp_mass(p.Tc, p)


def q_rx(r, p):
    return (-p.dH_rxn) * p.A_c * r


def q_cool(T, p):
    return p.perimeter * p.U * (T - p.Tc)


def dTdz(r, T, p):
    return (q_rx(r, p) - q_cool(T, p)) / heat_capacity_flow(p)
