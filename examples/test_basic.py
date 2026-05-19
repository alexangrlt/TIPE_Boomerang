from core import simulate_projectile, plot_trajectory_3d
from core import Boomerang_standard
import numpy as np

# Conditions initiales realistes d'un lancer de boomerang droitier :
# - 20 m/s en +X (~72 km/h)
# - Legere composante verticale +1 m/s
# - Hauteur de lancer 1.5 m
# - omega : 150 rad/s ~ 1430 rpm (valeur haute mais realiste pour un lancer fort)
# - Inclinaison Y -75 deg dans le simulateur (quasi-vertical)
result = simulate_projectile(
    np.array([0.0, 0.0, 1.5]),
    np.array([20.0, 0.0, 1.0]),
    Boomerang_standard,
    dt=0.0002,
    t_max=15,
)

# simulate_projectile retourne maintenant 6 valeurs
px, py, pz, pos, R_rot, omega_mag = result

plot_trajectory_3d(pos)

# Diagnostic : evolution du spin au cours du temps
import matplotlib.pyplot as plt
t = np.arange(len(omega_mag)) * 0.0002
plt.figure()
plt.plot(t, omega_mag)
plt.xlabel("Temps (s)")
plt.ylabel("||omega|| (rad/s)")
plt.title("Evolution du spin")
plt.grid(True)
plt.show()
