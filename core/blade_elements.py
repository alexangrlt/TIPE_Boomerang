import numpy as np
from core.boomerang_config import BoomerangConfig
from core.xflr_data import load_data


Cl_p1d, Cd_p1d = load_data("NACA_4420_T1_Re0.080_M0.00_N9.0.txt")

# Twist lineaire de la pale : 0 deg a la racine, -TWIST_TIP deg au bout.
#
# SIGNE NEGATIF : pour un boomerang droitier (rotation anti-horaire vue de dessus),
# le bout de pale avance vite -> alpha_aero y est deja grand.
# Un twist negatif au bout REDUIT alpha local pour eviter la sur-portance,
# ce qui equilibre la distribution de portance sur l'envergure.
# Cela correspond physiquement a un vrillage "nose-down" au bout (bord d'attaque
# qui descend), signe conventionnel negatif.
# Un twist positif ferait l'inverse : sur-portance au bout -> effondrement de la
# precession en fin de vol (trajectoire droite = bug observe).
TWIST_TIP = 14.0  # degres (valeur absolue, le signe est gere ci-dessous)


def get_blade_element(config):
    """
    Retourne les troncons de chaque pale avec twist lineaire negatif.

    Chaque troncon contient :
    - r          : position radiale (m)
    - dr         : largeur du troncon (m)
    - c          : corde locale (m)
    - dS         : surface du troncon (m^2)
    - angle_pale : angle de la pale dans le plan du boomerang (rad)
    - vect_unit  : vecteur unitaire de l'envergure (repere boomerang)
    - twist      : angle de vrillage local (rad), NEGATIF au bout (nose-down)

    Le desequilibre avancant/reculant emerge naturellement de
    v_rel = v_cm + cross(omega, r_vec) : la pale avancante a une vitesse
    relative plus grande => plus de portance => couple de precession.
    """
    r_values     = config.e
    dr           = r_values[1] - r_values[0]
    chord_values = np.linspace(config.c_root, config.c_tip, config.N_troncons)

    # Twist lineaire : 0 a la racine, -TWIST_TIP au bout (signe negatif = nose-down)
    twist_values = np.linspace(0.0, -np.radians(TWIST_TIP), config.N_troncons)

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
