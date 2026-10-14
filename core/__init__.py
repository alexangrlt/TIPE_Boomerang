"""Package principal de modélisation et simulation 3D du boomerang."""

from .basic_simulator import (
    plot_angles,
    plot_rot,
    plot_trajectory_3d,
    simulate_projectile,
)
from .blade_elements import POLAR_FILE, get_blade_element
from .boomerang_config import Boomerang_standard, BoomerangConfig
from .xflr_data import load_data

__all__ = [
    "BoomerangConfig",
    "Boomerang_standard",
    "POLAR_FILE",
    "get_blade_element",
    "load_data",
    "plot_angles",
    "plot_rot",
    "plot_trajectory_3d",
    "simulate_projectile",
]
