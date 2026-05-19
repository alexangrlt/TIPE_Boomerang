import numpy as np
from core.boomerang_config import BoomerangConfig
from core.xflr_data import load_data


Cl_p1d, Cd_p1d = load_data("NACA_4420_T1_Re0.080_M0.00_N9.0.txt")

# Twist lineaire de la pale : 0 deg a la racine, TWIST_TIP deg au bout.
# Sans twist, tous les elements ont le meme angle geometrique alors que
# la vitesse relative croit en r*omega -> sur-portance au bout, sous-portance
# a la racine. Le twist compense ca (meme principe qu'une helice).
TWIST_TIP = 12.0  # degres


def get_blade_element(config):
    """
    Retourne les troncons de chaque pale avec twist lineaire.

    Chaque troncon contient :
    - r          : position radiale (m)
    - dr         : largeur du troncon (m)
    - c          : corde locale (m)
    - dS         : surface du troncon (m^2)
    - angle_pale : angle de la pale dans le plan du boomerang (rad)
    - vect_unit  : vecteur unitaire de l'envergure (repere boomerang)
    - twist      : angle de vrillage local (rad), positif = bord d'attaque monte

    Le desequilibre avancant/reculant emerge naturellement de
    v_rel = v_cm + cross(omega, r_vec) : la pale avancante a une vitesse
    relative plus grande => plus de portance => couple de precession.
    """
    r_values     = config.e
    dr           = r_values[1] - r_values[0]
    chord_values = np.linspace(config.c_root, config.c_tip, config.N_troncons)

    # Twist lineaire : 0 a la racine, TWIST_TIP au bout
    twist_values = np.linspace(0.0, np.radians(TWIST_TIP), config.N_troncons)

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
                "twist":      twist_values[j],
            })
    return elements
