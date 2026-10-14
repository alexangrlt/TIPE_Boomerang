"""Simulation basique de la trajectoire 3D et de l'évolution du spin."""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from core import Boomerang_standard, plot_trajectory_3d, simulate_projectile

ASSETS_DIR = Path(__file__).resolve().parent.parent / "assets"


def main() -> None:
    position_init = np.array([0.0, 0.0, 1.5])
    vitesse_init = np.array([25.0, 0.0, -1.0])

    res = simulate_projectile(
        position_init,
        vitesse_init,
        config=Boomerang_standard,
        dt=0.0002,
        t_max=15.0,
    )

    # Affichage de la trajectoire 3D
    plot_trajectory_3d(res["pos"])

    # Affichage de la décroissance du spin
    t = np.arange(len(res["omega"])) * res["dt"]
    plt.figure()
    plt.plot(t, res["omega"])
    plt.xlabel("Temps (s)")
    plt.ylabel("||omega|| (rad/s)")
    plt.title("Evolution du spin")
    plt.grid(True)
    output_img = ASSETS_DIR / "spin_decay.png"
    output_img.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_img, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
