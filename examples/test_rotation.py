from scipy.spatial.transform import Rotation as R
import numpy as np


from core import simulate_projectile, plot_rot, plot_angles
from core import Boomerang_standard


px, py, pz, rotation = simulate_projectile([0, 0, 2], [5, 0, 2], Boomerang_standard)
plot_rot(rotation)

plot_angles(rotation, dt=0.01)