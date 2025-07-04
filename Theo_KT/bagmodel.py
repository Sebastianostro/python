# bagmodel.py

import numpy as np
from constants import hquerfm, B, Z0, alpha_c, x, m_light, m_strange
import scipy as sp

def energy_terms(R, n, ns, m=m_light, ms=m_strange, S=0, I=0):
    kinetic = n * np.sqrt(m**2 + (x / R)**2) + ns * np.sqrt(ms**2 + (x / R)**2)
    bag_energy = (B**4) * (4 / 3) * np.pi * R**3 - Z0 / R
    delta_em = (4 / 3) * alpha_c / R * (n * (n - 6) + S * (S + 1) + 3 * I * (I + 1))
    return kinetic + bag_energy 
#+ delta_em

def minimize_energy(n, ns, S, I, m=m_light, ms=m_strange):
    result = sp.optimize.minimize_scalar(lambda R: energy_terms(R, n, ns, m, ms, S, I), bounds=(0,10), method='bounded')
    R_opt = result.x
    M_min = result.fun
    return R_opt, M_min
