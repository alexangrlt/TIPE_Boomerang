"""Vérification et tracé des courbes polaires interpolées Cl(alpha) et Cd(alpha) du profil NACA 4420."""

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from core import POLAR_FILE, load_data


def main() -> None:
    # Chargement des interpolateurs depuis le fichier XFLR5
    Cl_func, Cd_func = load_data(POLAR_FILE)

    # Évaluation sur la plage d'incidence [-5°, 25°]
    alpha_test = np.linspace(-5.0, 25.0, 200)
    Cl_vals = Cl_func(alpha_test)
    Cd_vals = Cd_func(alpha_test)

    # Affichage côte à côte
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(alpha_test, Cl_vals, "b-")
    ax1.set_xlabel("Alpha (°)")
    ax1.set_ylabel("Cl")
    ax1.set_title("Cl(α) — NACA 4420")
    ax1.grid(True)

    ax2.plot(alpha_test, Cd_vals, "r-")
    ax2.set_xlabel("Alpha (°)")
    ax2.set_ylabel("Cd")
    ax2.set_title("Cd(α) — NACA 4420")
    ax2.grid(True)

    plt.tight_layout()

    output_img = Path(__file__).resolve().parent.parent / "assets" / "polar_check.png"
    output_img.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_img, dpi=200)
    plt.close(fig)

    # Vérification rapide sur quelques angles caractéristiques
    print(f"Cl à 5°  : {float(Cl_func(5)):.4f}")
    print(f"Cl à 10° : {float(Cl_func(10)):.4f}")
    print(f"Cd à 5°  : {float(Cd_func(5)):.4f}")


if __name__ == "__main__":
    main()
