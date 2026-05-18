import numpy as np
from core.boomerang_config import BoomerangConfig
from core.xflr_data import load_data


Cl_p1d, Cd_p1d = load_data("NACA_4420_T1_Re0.080_M0.00_N9.0.txt")

def get_blade_element(config):
    """
    Retourne les troncons de chaque pale.
    Chaque pale est definie par son angle (0 deg et 74 deg).
    Les troncons vont du centre (r=0.005) jusqu'au bout (r=R_pale).
    Le desequilibre avancant/reculant emerge naturellement de
    v_rel = v_cm + cross(omega, r_vec) : la pale qui avance
    (vitesse translation dans le meme sens que la rotation)
    a une vitesse relative plus grande => plus de portance.
    """
    r_values     = config.e
    dr           = r_values[1] - r_values[0]
    chord_values = np.linspace(config.c_root, config.c_tip, config.N_troncons)

    elements = []
    for angle_deg in config.angles_pales:
        angle_rad = np.radians(angle_deg)
        vect_unit = np.array([np.cos(angle_rad), np.sin(angle_rad), 0.0])
        for j, r in enumerate(r_values):
            c  = chord_values[j]
            dS = c * dr
            elements.append({
                "r":          r,
                "dr":         dr,
                "c":          c,
                "dS":         dS,
                "angle_pale": angle_rad,
                "vect_unit":  vect_unit,
            })
    return elements
