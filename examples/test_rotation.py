from scipy.spatial.transform import Rotation as R
import numpy as np


from core.basic_simulator import simulate_projectile, plot_rot


px, py, pz, rotation = simulate_projectile([0, 0, 2], [5, 0, 2])
plot_rot(rotation)