"""(7) Pressure drop - Phase 2 dimensionless Ergun parameterization (inertia-dominated form):

    dp/dzeta = -Pi (T/T_c)/p,   p = P/P_ref, zeta = z/L    <=>    dP/dz = -(Pi P_ref^2/L)(T/T_c)/P

Pi is the study parameter. No particle diameter, bed voidage or bulk density is used (none is
documented for the adopted dataset). This implements and verifies the MATHEMATICAL formulation;
the absolute industrial pressure drop is NOT independently validated.

Modes (Phase 2): 'none'     baseline: P = P_ref everywhere (benchmark)
                 'outlet'   P(L) = P_ref, inlet P0 from shooting (primary pressure-drop configuration)
                 'inlet'    P(0) = P_ref, only for Pi <= Pi_inlet_max (controlled comparison)
                 'constant' P = P0 everywhere (same inlet pressure, no decline; comparison)
Isothermal (T = T_c) relations: P(z)^2 = P0^2 - 2 Pi P_ref^2 z/L, hence
    outlet spec: P0 = P_ref sqrt(1 + 2 Pi);  inlet spec: P(L) = P_ref sqrt(1 - 2 Pi);
    Pi = [(P0/P_ref)^2 - (P(L)/P_ref)^2] / 2.
Non-isothermal identity: P(L)^2 = P0^2 - 2 Pi P_ref^2 int_0^1 (T/T_c) dzeta.
"""
import numpy as np

MODES = ("none", "outlet", "inlet", "constant")
PD_MODES = ("outlet", "inlet")


def coeff(p, Pi):
    return Pi * p.P_ref ** 2 / p.L


def dPdz(T, P, p, Pi):
    return -coeff(p, Pi) * (T / p.Tc) / P


def dPdz_partials(T, P, p, Pi):
    c = coeff(p, Pi)
    return -c / (p.Tc * P), c * (T / p.Tc) / P ** 2


def isothermal_profile(z, P0, p, Pi):
    return np.sqrt(P0 ** 2 - 2.0 * Pi * p.P_ref ** 2 * np.asarray(z) / p.L)


def isothermal_P0_outlet(p, Pi):
    return p.P_ref * np.sqrt(1.0 + 2.0 * Pi)


def isothermal_Pout_inlet(p, Pi):
    return p.P_ref * np.sqrt(1.0 - 2.0 * Pi)


def Pi_from_isothermal(P0, PL, p):
    return ((P0 / p.P_ref) ** 2 - (PL / p.P_ref) ** 2) / 2.0
