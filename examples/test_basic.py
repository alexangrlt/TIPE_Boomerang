from core import simulate_projectile, plot_trajectory_3d
from core import Boomerang_standard
import numpy as np

# Conditions initiales réalistes d'un lancer de boomerang
# Vitesse de lancer ~22 m/s vers l'avant, légère composante verticale
# Hauteur de lancer ~1.5 m

px, py, pz, pos, R_rot = simulate_projectile(
    np.array([0.0, 0.0, 1.5]),   # position_init : hauteur de lancer réaliste
    np.array([22.0, 0.0, 2.0]),  # vitesse_init : ~22 m/s vers l'avant
    Boomerang_standard,
    dt=0.0002,
    t_max=15,
)

plot_trajectory_3d(pos)
