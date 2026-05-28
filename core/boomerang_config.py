"""
Configuration physique du boomerang
"""

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)  # immutables
class BoomerangConfig:
    """Settings Boomerang 2 pales"""

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

    @property
    def envergure(self):
        """envergure du boomerang"""
        return self.L * np.sqrt(
            2
        )  # sqrt(2) car on considère un carré et la diag d'un carré est la longueur *sqrt(2), en m

    def surface(self):
        """surface totale du boomerang"""
        return (self.c_root + self.c_tip) * self.R_pale  # calcul de surface, en m²

    def matrice_inertie(self):
        """
        Tenseur d'inertie du boomerang : 2 pales rectangulaires.

        CORRECTION (Steiner) : le tenseur de chaque pale est d'abord calcule dans
        le repere de la pale (origine = centre de masse de la pale, non du boomerang),
        puis transpose au centre du boomerang via le theoreme d'Huygens-Steiner :
            I_total = I_cm_pale + m * (|d|^2 * Id - d dyadic d)
        ou d est le vecteur entre le centre de masse de la pale et le centre du boomerang.
        """
        m_pale = self.masse / 2          # masse d'une pale
        L = self.R_pale                  # longueur de pale = R_pale (extension radiale)
        c = (self.c_root + self.c_tip) / 2  # corde moyenne

        # Inertie dans le repere propre de la pale (axe pale = axe Y local)
        # Pale rectangulaire : Ixx = 1/12*m*c^2, Iyy = 1/3*m*L^2 (Huygens par rapport a la racine)
        # Mais ici on veut I par rapport au CM de la pale, situe a L/2 de la racine :
        # Ixx_cm = 1/12*m*c^2  (invariant : la corde est centree)
        # Iyy_cm = 1/12*m*L^2  (barre mince : I_cm = 1/12*m*L^2)
        # Izz_cm = Ixx_cm + Iyy_cm (theoreme de König pour profil mince)
        Ixx_cm = 1.0/12.0 * m_pale * c**2
        Iyy_cm = 1.0/12.0 * m_pale * L**2
        Izz_cm = Ixx_cm + Iyy_cm
        I_pale_cm = np.diag([Ixx_cm, Iyy_cm, Izz_cm])

        I_total = np.zeros((3, 3))

        for angle_deg in self.angles_pales:
            angle_rad = np.radians(angle_deg)

            # Matrice de rotation autour de z (de la pale dans le repere boomerang)
            Rz = np.array([
                [ np.cos(angle_rad), -np.sin(angle_rad), 0],
                [ np.sin(angle_rad),  np.cos(angle_rad), 0],
                [ 0,                  0,                  1]
            ])

            # Rotation du tenseur d'inertie au CM de la pale vers le repere boomerang
            I_pale_rot = Rz @ I_pale_cm @ Rz.T

            # FIX bug 8 : le CM de la pale est a R_pale/2 de l'origine, pas L/2.
            # L (0.195 m) est la longueur geometrique totale du boomerang deplie,
            # tandis que R_pale (0.112 m) est l'extension radiale effective de chaque pale.
            # Utiliser L surestimait le bras de levier de Steiner d'un facteur ~1.7
            # et donc le moment d'inertie, ce qui amortissait trop la precession.
            d = Rz @ np.array([self.R_pale / 2.0, 0.0, 0.0])

            # Theoreme de Steiner : I_origine = I_cm + m*(|d|^2*Id3 - d x d^T)
            d_sq = np.dot(d, d)
            steiner = m_pale * (d_sq * np.eye(3) - np.outer(d, d))

            I_total += I_pale_rot + steiner

        return I_total

Boomerang_standard = BoomerangConfig()
