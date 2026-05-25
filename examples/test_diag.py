"""Script de diagnostic - affiche vitesse, forces et inclinaison au cours du temps.
Permet d'identifier pourquoi le boomerang ne revient pas (x fixe en fin de vol).
"""
from core import simulate_projectile
from core import Boomerang_standard
import numpy as np
import matplotlib.pyplot as plt

res = simulate_projectile(
    np.array([0.0, 0.0, 1.5]),
    np.array([20.0, 0.0, 1.0]),
    Boomerang_standard,
    dt=0.0002,
    t_max=15,
)

pos     = np.array(res["pos"])
vit     = res["vitesse"]
F_aero  = res["F_aero"]
omega   = res["omega"]
incl    = res["incl"]
dt      = res["dt"]
t       = np.arange(len(pos)) * dt

# --- Vue de dessus (X-Y) ---
fig, axes = plt.subplots(2, 3, figsize=(15, 8))
fig.suptitle("Diagnostic trajectoire boomerang", fontsize=13)

axes[0,0].plot(pos[:,0], pos[:,1])
axes[0,0].plot(pos[0,0], pos[0,1], 'go', label='depart')
axes[0,0].plot(pos[-1,0], pos[-1,1], 'ro', label='arrivee')
axes[0,0].set_xlabel("X (m)"); axes[0,0].set_ylabel("Y (m)")
axes[0,0].set_title("Trajectoire vue de dessus (X-Y)")
axes[0,0].legend(); axes[0,0].grid(True); axes[0,0].set_aspect('equal')

# --- Vitesses ---
axes[0,1].plot(t, vit[:,0], label='vx')
axes[0,1].plot(t, vit[:,1], label='vy')
axes[0,1].plot(t, vit[:,2], label='vz')
axes[0,1].set_xlabel("Temps (s)"); axes[0,1].set_ylabel("v (m/s)")
axes[0,1].set_title("Composantes de vitesse")
axes[0,1].legend(); axes[0,1].grid(True)

# --- Norme vitesse horizontale ---
V_horiz = np.sqrt(vit[:,0]**2 + vit[:,1]**2)
axes[0,2].plot(t, V_horiz)
axes[0,2].set_xlabel("Temps (s)"); axes[0,2].set_ylabel("|v_horiz| (m/s)")
axes[0,2].set_title("Vitesse horizontale (freinée si trop forte trainee)")
axes[0,2].grid(True)

# --- Forces aerodynamiques ---
axes[1,0].plot(t, F_aero[:,0], label='Fx_aero')
axes[1,0].plot(t, F_aero[:,1], label='Fy_aero')
axes[1,0].plot(t, F_aero[:,2], label='Fz_aero')
axes[1,0].set_xlabel("Temps (s)"); axes[1,0].set_ylabel("F (N)")
axes[1,0].set_title("Forces aerodynamiques (repere monde)")
axes[1,0].legend(); axes[1,0].grid(True)

# --- Inclinaison du plan ---
axes[1,1].plot(t, incl)
axes[1,1].axhline(90, color='r', linestyle='--', label='plan horizontal (90 deg)')
axes[1,1].set_xlabel("Temps (s)"); axes[1,1].set_ylabel("inclinaison (deg)")
axes[1,1].set_title("Inclinaison du plan / vertical (0=vertical, 90=horizontal)")
axes[1,1].legend(); axes[1,1].grid(True)

# --- Spin ---
axes[1,2].plot(t, omega)
axes[1,2].set_xlabel("Temps (s)"); axes[1,2].set_ylabel("omega (rad/s)")
axes[1,2].set_title("Evolution du spin")
axes[1,2].grid(True)

plt.tight_layout()
plt.show()

# Impression des valeurs en fin de trajectoire
print(f"--- Fin de trajectoire ---")
print(f"Position finale  : x={pos[-1,0]:.2f} m  y={pos[-1,1]:.2f} m  z={pos[-1,2]:.2f} m")
print(f"Vitesse finale   : vx={vit[-1,0]:.2f}  vy={vit[-1,1]:.2f}  vz={vit[-1,2]:.2f} m/s")
print(f"Inclinaison finale (0=vert, 90=horiz) : {incl[-1]:.1f} deg")
print(f"Spin final       : {omega[-1]:.1f} rad/s")
print(f"Duree de vol     : {t[-1]:.2f} s")

# Identification du moment ou vx devient proche de 0
idx_vx0 = np.argmin(np.abs(vit[:,0]))
print(f"\nvx ~ 0 a t={t[idx_vx0]:.2f}s  (position x={pos[idx_vx0,0]:.2f} m  incl={incl[idx_vx0]:.1f} deg)")
