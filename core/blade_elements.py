import numpy as np
from core.boomerang_config import BoomerangConfig
from core.xflr_data import load_data


Cl_p1d, Cd_p1d = load_data("NACA_4420_T1_Re0.080_M0.00_N9.0.txt")

def get_blade_element(config):
    r_values = config.e  # init des rayons avec le linspace
    dr = r_values[1] - r_values[0]  # pas uniforme

    chord_values = np.linspace(
        config.c_root, config.c_tip, config.N_troncons
    )  # corde pour chaque r

    elements = []
    for i in config.angles_pales:
        angle_rad = np.radians(i)
        # vecteur unitaire dans la direction de la pale
        vect_unit = np.array([np.cos(angle_rad), np.sin(angle_rad), 0])
        for j, r in enumerate(r_values):
            c = chord_values[j]
            dS = c * dr
            base = {
                "dr": dr,
                "c": c,
                "dS": dS,
                "angle_pale": angle_rad,
                "vect_unit": vect_unit,
            }
            # tronçon coté + (bout de pale)
            elements.append({**base, "r": r})
            # tronçon coté - (base de pale de l'autre côté du centre)
            # nécessaire pour le déséquilibre avançant/reculant
            elements.append({**base, "r": -r})
    return elements
