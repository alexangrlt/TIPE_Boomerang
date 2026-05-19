import numpy as np
from scipy.interpolate import interp1d


def load_data(chemin_d_acces_au_fichier):
    """
    Charge les donnees polaires XFLR5 et retourne des interpolateurs Cl(alpha)
    et Cd(alpha).

    Attention : les fichiers XFLR5 ont parfois des trous dans la plage d'alpha
    (ex: donnees absentes entre 3 et 9.5 deg). On remplit ces trous par
    interpolation lineaire avant de construire l'interpolateur final, pour
    eviter un saut brutal de Cl qui fausserait le calcul de portance.
    """
    data = np.loadtxt(chemin_d_acces_au_fichier, skiprows=11)

    # Recup alpha, Cl, Cd
    alpha = data[:, 0]
    Cl    = data[:, 1]
    Cd    = data[:, 2]

    # Tri par alpha croissant (securite)
    ordre = np.argsort(alpha)
    alpha = alpha[ordre]
    Cl    = Cl[ordre]
    Cd    = Cd[ordre]

    # Remplissage des trous : on cree une grille reguliere de -5 a 25 deg
    # par pas de 0.5 deg et on interpole lineairement sur les points existants
    # avant de construire l'interpolateur final.
    alpha_dense = np.arange(alpha[0], alpha[-1] + 0.1, 0.5)
    Cl_dense = np.interp(alpha_dense, alpha, Cl)
    Cd_dense = np.interp(alpha_dense, alpha, Cd)

    # Interpolateurs finaux avec extrapolation constante hors plage
    Cl_p1d = interp1d(alpha_dense, Cl_dense,
                      kind='linear', bounds_error=False,
                      fill_value=(Cl_dense[0], Cl_dense[-1]))
    Cd_p1d = interp1d(alpha_dense, Cd_dense,
                      kind='linear', bounds_error=False,
                      fill_value=(Cd_dense[0], Cd_dense[-1]))

    return Cl_p1d, Cd_p1d
