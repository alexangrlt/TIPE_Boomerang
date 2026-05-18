import matplotlib.pyplot as plt
import numpy as np
from core.xflr_data import load_data

# Charge les données
Cl_func, Cd_func = load_data("NACA_4420_T1_Re0.080_M0.00_N9.0.txt")

# Teste sur une plage d'alpha
alpha_test = np.linspace(-5, 25, 200)
Cl_vals = Cl_func(alpha_test)
Cd_vals = Cd_func(alpha_test)

# Affiche
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(alpha_test, Cl_vals, 'b-')
ax1.set_xlabel("Alpha (°)")
ax1.set_ylabel("Cl")
ax1.set_title("Cl(α) — NACA 4420")
ax1.grid(True)

ax2.plot(alpha_test, Cd_vals, 'r-')
ax2.set_xlabel("Alpha (°)")
ax2.set_ylabel("Cd")
ax2.set_title("Cd(α) — NACA 4420")
ax2.grid(True)

plt.tight_layout()
plt.savefig("polar_check.png")
plt.show()

# Vérif rapide sur quelques valeurs
print(f"Cl à 5° : {Cl_func(5):.4f}")
print(f"Cl à 10° : {Cl_func(10):.4f}")
print(f"Cd à 5° : {Cd_func(5):.4f}")