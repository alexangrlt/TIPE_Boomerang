from scipy.spatial.transform import Rotation as R
import numpy as np


from core import simulate_projectile, plot_rot, plot_angles
from core import Boomerang_standard


px, py, pz, pos, rotation = simulate_projectile(
    [0, 0, 2], [10, 0, 5], Boomerang_standard
)
plot_rot(rotation)

plot_angles(rotation, dt=0.01)
