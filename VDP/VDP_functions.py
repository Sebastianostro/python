# This file contains various functions related to particle physics calculations, particularly focusing on the Bethe-Bloch formula and related concepts.
import math
# Constants
import VDP_constants as const
# Coulomb correction function for electron/positron
def f_coul(Z):
    a = Z/137
    f = a**2*((1+a**2)**(-1)+0.20206-0.0369*a**2+0.0083*a**4-0.002*a**6)
    return f

# Separate helper function for the common factor of the Bethe-Bloch formula
def bethe_bloch_common_factor():
    # Calculate the common factor
    K = 2 * math.pi * const.N_A * const.r_e**2 * const.m_e
    return K

# Bethe-Bloch formula for the energy loss
def bethe_bloch(Z_det, A_det, rho_det, Z_projectile, beta_projectile, W_max, I_potential, density_correction, coulomb_correction):
    # Calculate beta_gamma
    beta_gamma = beta_projectile / math.sqrt(1 - beta_projectile**2)
    # Calculate common factor
    K = 2 * math.pi * const.N_A * const.r_e**2 * const.m_e
    dE = K * Z_det / A_det * rho_det * Z_projectile**2 / beta_projectile**2 * (math.log(2 * const.m_e * beta_gamma**2 * W_max / I_potential**2) - 2 * beta_projectile**2 - density_correction - 2 * coulomb_correction/Z_det)
    return dE

# number of atoms in 1 cm^3
def rho_to_N(rho, A):
    # Calculate number of atoms per cm^3
    N = rho * const.N_A / A
    return N
# Sollte 2,7*10^19 sein

# Function for the radiative correction
def phi_rad(Z, E0, shield_ind):
    # Calculate common factor
    cf = 4 * Z**2 * const.r_e**2 * const.alpha
    # Calculate the radiative correction
    if shield_ind == 0:
        phi = cf * (math.log(2 * E0 / const.m_e) - 1/3 - f_coul(Z))
    elif shield_ind == 1:
        phi = cf * (math.log(183*Z**(-1/3)) + 1/18 - f_coul(Z))
    return phi

# Function to check for the degree of shielding
def check_phi_rad(Z):
    # Check if the radiative correction is applicable
    limit = const.alpha * const.m_e * Z**(1/3)
    return limit

# Function for the energy loss
def dE_dx(E0, Z, A, rho, shield_ind):
    # Calculate the energy loss
    dE = rho_to_N(rho, A) * E0 * phi_rad(Z, E0, shield_ind)
    return dE

# Approximation of the critical energy for electrons
def E_crit_approx(Z):
    # Calculate the critical energy using an approximation
    E_crit = 800/(Z+1.2)
    return E_crit

# Simplified radiation length function
def L_rad(Z, A):
    # Calculate the radiation length
    X0 = 716.4 * A / (Z*(Z+1) * math.log(287/math.sqrt(Z)))
    return X0

# Bethe-Heitler distribution function
def bh_distribution(E, E0, Z, A, material_length):
    # Calculate the length in terms of the radiation length
    t = material_length / L_rad(Z, A)
    # Calculate the relative energy
    x = E / E0
    # Calculate the distribution
    dN = (-math.log(x))**(t/math.log(2)-1)/math.gamma(t/math.log(2))
    return dN 

# Function to calculate the maximum energy loss (or energy transfer)
def E_max(beta_gamma, gamma, M):
    # Calculate the maximum energy transfer under the assumption of natural units
    # where c = 1 and masses are in MeV/c^2
    E_max = 2 * const.m_e * beta_gamma**2 / (1 + 2 * gamma * (const.m_e / M) + (const.m_e / M)**2)
    return E_max