import numpy as np
from core.boomerang_config import BoomerangConfig

def get_blade_element(config):
    r_values=config.e               #init des rayons avec le linspace
    dr = r_values[1]-r_values[0]    #le pas de 2 rayons consécutifs -- uniforme
    
    chord_values=np.linspace(config.c_root,config.c_tip, config.N_troncons)     #valeurs des cordes pour chaque r
    
    elements=[]
    for i in config.angles_pales:
        angle_rad=np.rad(i)
        """je me crée un vecteur unitaire associé à cet angle"""
        vect_unit=np.array([np.cos(angle_rad),np.sin(angle_rad),0]) #[x,y,z]
        for j,r in enumerate(r_values):     #enumerate donne l'indice ET la valeur !
            c=chord_values[j]
            dS=c*dr
            elements.append({"r": r,
                "dr": dr,
                "chord": c,
                "dS": dS,
                "angle_pale": angle_rad,
                "u_rad": vect_unit,
            })
    return elements