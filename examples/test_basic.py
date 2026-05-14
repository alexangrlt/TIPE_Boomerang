from core import simulate_projectile, plot_trajectory_3d
from core import Boomerang_standard

# je suppose pour l'instant que le boomerang est lancé a partir de 2m

px, py, pz, pos, R_rot = simulate_projectile(
    [0, 0, 2],  # position_init
    [8, 0, 1],  # vitesse_init
    Boomerang_standard,
    dt=0.001,
    t_max=10,
)

plot_trajectory_3d(pos)
