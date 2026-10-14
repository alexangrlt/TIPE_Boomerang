"""Visualisation de l'évolution de l'orientation 3D et des angles d'Euler (roulis, tangage, lacet)."""

from core import Boomerang_standard, plot_angles, plot_rot, simulate_projectile


def main() -> None:
    res = simulate_projectile(
        [0.0, 0.0, 2.0],
        [10.0, 0.0, 5.0],
        config=Boomerang_standard,
    )

    rotation = res["rot"]
    dt = res["dt"]

    plot_rot(rotation)
    plot_angles(rotation, dt=dt)


if __name__ == "__main__":
    main()
