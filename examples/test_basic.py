from core.basic_simulator import simulate_projectile, plot_trajectory_3d


#je suppose pour l'instant que le boomerang est lancé a partir de 2m
#et une v_init de 5m/s selon x, 2m/s selon z (lancer parfait)

px, py, pz, R_rot = simulate_projectile([0, 0, 2], [5, 0, 2])
plot_trajectory_3d(px, py, pz,)