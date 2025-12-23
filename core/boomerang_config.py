"""
Configuration physique du boomerang
"""

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)  # immutables
class BoomerangConfig:
    """Settings Boomerang 2 pales (rectangulaires)"""

    # Parametres Geometriques
    L = 0.195  # longueur de pale en m
    l = 0.034  # largeur de pale en m
    h = 0.0031  # epaisseur de pale en m

    # Masse
    masse = 0.028  # masse en kg #masse de mon 4eme boomerang (1000%, 100% infill (triangles), 3h40 d'impression, PLA)

    # Parametres Aerodynamiques
    rho = 1.225  # densité de l'air en kg/m³
    Cx = 1.5  # coefficient de traînée                         #Je ne trouve aucune source donnant Cx (pour tester) donc en premiere approche j'ai pris le Cx de Corentin
    Cz = 0.45  # coefficient de portance initial               #https://www.math.uci.edu/~eesser/papers/justboom.pdf  ##y a un calcul a faire (ou plutot une experience....)
    Cm = 0.05  # coefficient de Magnus          #0.05 juste pr lancer la simu mais :je n'ai trouvé nulle part une valeur qui pourrait correspondre a mon boomerang et/ou un calcul pour arriver a cette valeur car elle se trouve par l'experience...

    @property
    def envergure(self):
        """envergure du boomerang"""
        return self.L * np.sqrt(
            2
        )  # sqrt(2) car on considère un carré et la diag d'un carré est la longueur *sqrt(2), en m

    def surface(self):
        """surface totale du boomerang"""
        return self.L * self.l * 2  # calcul de surface, en m²


Boomerang_standard = BoomerangConfig()
