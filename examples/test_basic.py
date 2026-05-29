from core import simulate_projectile, plot_trajectory_3d
from core import Boomerang_standard
import numpy as np
import matplotlib.pyplot as plt

res = simulate_projectile(
    np.array([0.0, 0.0, 1.5]),
    np.array([25.0, 0.0, 5.0]),
    Boomerang_standard,
    dt=0.0001,
    t_max=15,
)

plot_trajectory_3d(res["pos"])

t = np.arange(len(res["omega"])) * res["dt"]
plt.figure()
plt.plot(t, res["omega"])
plt.xlabel("Temps (s)")
plt.ylabel("||omega|| (rad/s)")
plt.title("Evolution du spin")
plt.grid(True)
plt.show()
