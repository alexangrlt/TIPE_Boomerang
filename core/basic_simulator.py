import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import (
    interp1d,
)  # je vais en avoir besoin pr remplir le vide entre les valeurs que va me donner XFLR
from core.boomerang_config import Boomerang_standard, BoomerangConfig

# ! mettre en place les tronçons, appliquer les forces dessus + calculs moments
# mais je galere vrmt a decouper le boomerang en troncon...
# prcq il n'est pas plein dans le sens où si je fais des petits cercles depuis le centre, j'ai pas mal de vide a des moments...
# ? mais mon idée à l'air d'être un "blade element theory" https://en.wikipedia.org/wiki/Blade_element_theory
# je vais voir si j'ai une idée sinon je demande de l'aide
# je me suis précipité sur cette idée mais je bloque...


def simulate_projectile(position_init, vitesse_init, config, dt=0.001, t_max=20):
    position = position_init.copy()
    vitesse = vitesse_init.copy()
    t = 0
    Px, Py, Pz = [], [], []
    pos = []
    R_rot, R_omega = [], []
    g = 9.81

    # ? j'ai considéré le boomerang a 25° p/r a la normale au sol (z)
    # ? donc 65° p/r a x (avec x vers l'avant, y vers la gauche et z vers le haut)

    rot_current = R.from_rotvec([65 * np.pi / 180, 0, 0])
    omega = np.array([0, 8, 0])  # on considère une vitesse angulaire constante ici

    while t < t_max and position[2] > 0:
        pos.append(
            [position[0], position[1], position[2]]
        )  # liste de listes des positions
        Px.append(position[0])
        Py.append(position[1])
        Pz.append(position[2])
        """Forces"""
        F_gravite = np.array([0, 0, -config.masse * g])
        # F_magnus = config.Cm * np.cross(omega, vitesse)
        # F_portance = (
        #     0.5
        #     * config.rho
        #     * np.linalg.norm(vitesse) ** 2
        #     * Boomerang_standard.surface()
        #     * config.Cz
        # ) * rot_current.apply(
        #     [0, 0, 1]
        # )  # prcq ça dépend de l'inclinaison (donc rotation) du boomerang

        # if np.linalg.norm(vitesse) > 0:
        #     F_trainee = (
        #         -0.5
        #         * config.rho
        #         * np.linalg.norm(vitesse) ** 2
        #         * Boomerang_standard.surface()
        #         * config.Cx
        #         * (
        #             vitesse / np.linalg.norm(vitesse)
        #         )  # pour avoir le sens de la force (le vecteur)
        #     )

        # else:
        #     F_trainee = np.array([0, 0, 0])

        # F_tot = F_gravite + F_magnus + F_portance + F_trainee
        F_tot = F_gravite  #! +  F_totPale
        # ? je dois faire un F_totPale qui est la somme de mes F_troncon.
        # ? sauf que je galere a faire mon decoupage en troncons...
        acceleration = F_tot / config.masse

        vitesse += acceleration * dt

        position = [
            position[0] + vitesse[0] * dt,
            position[1] + vitesse[1] * dt,
            position[2] + vitesse[2] * dt,
        ]

        R_rot.append(rot_current.as_rotvec())
        rot_increment = R.from_rotvec(omega * dt)
        rot_current = rot_current * rot_increment
        t += dt

    rotation = np.array(R_rot)

    # print(F_portance, F_trainee)
    print(pos)

    return Px, Py, Pz, pos, rotation


def plot_rot(rotation, title="Position Angulaire"):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot(rotation[:, 0], rotation[:, 1], rotation[:, 2], label="rotation")
    ax.set_xlabel("X (rad)")
    ax.set_ylabel("Y (rad)")
    ax.set_zlabel("Z (rad)")
    ax.set_title(title)
    ax.legend()
    plt.show()


def plot_trajectory_3d(pos, title="Trajectoire du Projectile"):
    pos = np.array(pos)
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot(pos[:, 0], pos[:, 1], pos[:, 2])
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_zlabel("Z (m)")
    ax.set_title(title)
    plt.show()


def plot_angles(rotation, dt):
    # Conversion de la liste rotation en 'objet' rotations
    rotations = [R.from_rotvec(vec) for vec in rotation]

    # Conversion de l'objet rotations en angle dans xyz (orientation)
    angles = np.array([r.as_euler("xyz", degrees=True) for r in rotations])
    nb_points = len(rotation)
    t = np.arange(nb_points) * dt

    plt.figure()
    plt.plot(t, angles[:, 0], label="Angle X (roulis)")
    plt.plot(t, angles[:, 1], label="Angle Y (tangage)")
    plt.plot(t, angles[:, 2], label="Angle Z (lacet)")

    plt.xlabel("Temps (en s)")
    plt.ylabel("Angle (en degrés)")
    plt.legend()
    plt.grid(True)
    plt.show()
