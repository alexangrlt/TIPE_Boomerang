"""
Configuration physique du boomerang
"""

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)  # immutables
class BoomerangConfig:
    """Settings Boomerang 2 pales (rectangulaires)"""

    # Parametres Geometriques
    L = 0.195  # longueur de pale en m                                                                      #? mesuré
    l = 0.034  # largeur de pale en m                                                                       #? mesuré
    h = 0.0031  # epaisseur de pale en m                                                                    #? mesuré

    R_pale = 0.112  # rayon (centre boomerang jusqu'aux extremités des pales) en m                          #? mesuré
    rho_boomerang = 1.24  # densité du filament en g/cm³ donc en kg/dm³ donc kg/L                           #? donnée vendeur (cura)
    N_troncons = 20  # ? arbitraire
    angles_pales = (
        0,
        74,
    )  # angle de la 1ere pale et de la deuxieme en degré                                                   #? par calcul avec trigos et dessin du boomerang
    c_root = 0.042  # largeur de la racine (à 0.0015 m du bout) en m #chord                                 #? mesuré
    c_tip = 0.033  # largeur du bout des pales (à 0.001 m du bout) en m #chord                              #? mesuré

    e = np.linspace(
        0.005,  # on enleve l'epaisseur de debut, d'ou le 0.005 (sinon on va avoir une epaisseur infinie)
        R_pale,
        N_troncons,  # nb de valeurs
    )

    # Masse
    masse = 0.028  # masse en kg #4eme boomerang (1000%, 100% infill (triangles), 3h40 d'impression, PLA)   #? donnée par le slicer u(m)=0.002        #faudrait que je le pèse réellement mais pour le code je pense pas que ça change grand chose...

    # Parametres Aerodynamiques
    rho_air = 1.225  # densité de l'air en kg/m³
    # Cx = 1.5  # coefficient de traînée                         #Je ne trouve aucune source donnant Cx (pour tester) donc en premiere approche j'ai pris le Cx de Corentin
    # Cz = 0.45  # coefficient de portance initial               #https://www.math.uci.edu/~eesser/papers/justboom.pdf  ##y a un calcul a faire (ou plutot une experience....)
    # Cm = 0.05  # coefficient de Magnus          #0.05 juste pr lancer la simu mais :je n'ai trouvé nulle part une valeur qui pourrait correspondre a mon boomerang et/ou un calcul pour arriver a cette valeur car elle se trouve par l'experience...

    @property
    def envergure(self):
        """envergure du boomerang"""
        return self.L * np.sqrt(
            2
        )  # sqrt(2) car on considère un carré et la diag d'un carré est la longueur *sqrt(2), en m

    def surface(self):
        """surface totale du boomerang"""
        return (self.c_root + self.c_tip) * self.R_pale  # calcul de surface, en m²


Boomerang_standard = BoomerangConfig()
