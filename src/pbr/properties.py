"""(2) Thermophysical/property functions.

Ideal gas. Total molar flow is constant (Delta n = 0), so the mixture molar mass stays M_wf.
c_p is constant (Phase 2 assumption). No property data beyond the parameter file are used.
"""


def molar_concentration(T, P, p):
    """Total molar concentration P/(R T) [mol/m3]."""
    return P / (p.R_gas * T)


def gas_density(T, P, p):
    """P M_wf/(R T) [kg/m3]."""
    return P * p.M_wf / (p.R_gas * T)


def volumetric_flow(T, P, p):
    """N_f R T / P [m3/s] per tube."""
    return p.N_f * p.R_gas * T / P


def superficial_velocity(T, P, p):
    """Q / A_c [m/s]."""
    return volumetric_flow(T, P, p) / p.A_c


def cp_mass(T, p):
    """Mass heat capacity [J/(kg K)]; constant in the Phase 2 model."""
    return p.cp
