"""Chargement et interpolation des polaires aérodynamiques issues de XFLR5."""

from pathlib import Path
import numpy as np
from scipy.interpolate import interp1d


def charger_polaire_xflr5(chemin_fichier: str | Path):
    """Charge les données polaires XFLR5 (format texte) et retourne

    des interpolateurs Cl(alpha) et Cd(alpha) sécurisés.

    Les angles d'incidence du fichier (en degrés) sont convertis
    en radians pour homogénéiser avec le reste du code physique.
    """
    path = Path(chemin_fichier)
    if not path.is_file():
        raise FileNotFoundError(f"Fichier de polaire introuvable : {path}")

    # Chargement en ignorant les lignes d'en-tête (souvent 11 lignes chez XFLR5)
    data = np.loadtxt(path, skiprows=11)

    # Extraction des colonnes (alpha en deg, Cl, Cd)
    alpha_deg = data[:, 0]
    Cl = data[:, 1]
    Cd = data[:, 2]

    # Sécurité : tri par angle d'incidence croissant
    ordre = np.argsort(alpha_deg)
    alpha_deg = alpha_deg[ordre]
    Cl = Cl[ordre]
    Cd = Cd[ordre]

    # Remplissage des trous éventuels via une grille dense en degrés
    alpha_dense_deg = np.arange(alpha_deg[0], alpha_deg[-1] + 0.1, 0.5)
    Cl_dense = np.interp(alpha_dense_deg, alpha_deg, Cl)
    Cd_dense = np.interp(alpha_dense_deg, alpha_deg, Cd)

    # Conversion de la grille d'incidence en radians pour le solveur
    alpha_dense_rad = np.deg2rad(alpha_dense_deg)

    # Construction des interpolateurs finaux (en radians)
    Cl_interpolateur = interp1d(
        alpha_dense_rad,
        Cl_dense,
        kind="linear",
        bounds_error=False,
        fill_value=(Cl_dense[0], Cl_dense[-1]),
    )
    Cd_interpolateur = interp1d(
        alpha_dense_rad,
        Cd_dense,
        kind="linear",
        bounds_error=False,
        fill_value=(Cd_dense[0], Cd_dense[-1]),
    )

    return Cl_interpolateur, Cd_interpolateur


def load_data(chemin_fichier: str | Path):
    """Alias anglais conservé pour les appelants existants."""
    return charger_polaire_xflr5(chemin_fichier)
