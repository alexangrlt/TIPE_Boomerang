from core import simulate_projectile, plot_trajectory_3d
from core import Boomerang_standard
import numpy as np

# Conditions initiales realistes d'un lancer de boomerang droitier :
# - Vitesse de lancer ~20 m/s vers l'avant (+X)  => ~72 km/h, typique
# - Legere composante verticale positive (~2 m/s) pour lancer legerement vers le haut
# - Position initiale : 1.5 m de hauteur (hauteur des mains)
# - omega initial : ~125 rad/s ~ 1200 rpm (tourne dans le sens +Z corps)
#   apres rotation initiale, exprime dans le monde par simulate_projectile
#
# L'inclinaison initiale (80 deg par rapport sol, rotation Y dans simulator)
# est coherente avec le lancer en +X : le plan du boomerang est quasi vertical
# et legrement incline vers la gauche (favorise la precession vers la gauche)
px, py, pz, pos, R_rot = simulate_projectile(
    np.array([0.0, 0.0, 1.5]),    # position_init : hauteur de lancer
    np.array([20.0, 0.0, 2.0]),   # vitesse_init  : lancer droit, legerement vers le haut
    Boomerang_standard,
    dt=0.0002,
    t_max=15,
)

plot_trajectory_3d(pos)
