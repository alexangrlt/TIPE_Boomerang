from core import simulate_projectile, plot_trajectory_3d
from core import Boomerang_standard
import numpy as np

# Conditions initiales realistes d'un lancer de boomerang droitier
# Vitesse de lancer ~15 m/s vers l'avant (+X)
# Pas de composante laterale artificielle
px, py, pz, pos, R_rot = simulate_projectile(
    np.array([0.0, 0.0, 1.5]),   # position_init : hauteur de lancer
    np.array([15.0, 0.0, 1.0]),  # vitesse_init : lancer droit
    Boomerang_standard,
    dt=0.0002,
    t_max=15,
)

plot_trajectory_3d(pos)
