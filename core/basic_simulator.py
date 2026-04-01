import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import (
    interp1d,
)  # je vais en avoir besoin pr remplir le vide entre les valeurs que va me donner XFLR
from core.boomerang_config import Boomerang_standard, BoomerangConfig
from core.blade_elements import get_blade_element

# ! mettre en place un compute_moments_be
# ? en faisant cela, je recup mon M_tot (moment total) de chaque element de pale et je l'intègre de la meme manière que les Forces (F_tot)


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

    elements = get_blade_element(config)

    while t < t_max and position[2] > 0:
        pos.append(
            [position[0], position[1], position[2]]
        )  # liste de listes des positions
        Px.append(position[0])
        Py.append(position[1])
        Pz.append(position[2])
        """Forces"""
        F_gravite = np.array([0, 0, -config.masse * g])

        F_pale = compute_forces_be(elements, vitesse, omega, rot_current, config)

        F_tot = F_gravite + F_pale

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
    # print(pos)

    return Px, Py, Pz, pos, rotation


def compute_forces_be(elements, v_translation, omega, rot_current, config):
    """calcule la force totale F_tot exercée sur chaque éléments

    Args:
        elements (dict): dictionnaireayant toutes les infos de chaque element du boomerang
        v_translation (array): tableau des vitesses de translation du centre de masse du boomerang en fonction du temps en 3axes
        omega (array): tableau des vitesses de rotation 3axes en fonction du temps
        rot_current (scipy rot): orientation du boomerang
        config (_type_): appel des données de config
    """
    Cx_temp = 0.8  # ? a changer avec les données de xflr
    Cz_temp = 0.45  # ? a changer avec les données de xflr
    F_tot = np.zeros(
        3
    )  # init ma liste des forces (3axes) avec des zeros, je change apres les valeurs
    for e in elements:
        # position du troncon dans le repère terrestre (absolu)
        "Toutes les infos dans le ref absolu je les noterai ..._abs"
        vect_unit_abs = rot_current.apply(e["vect_unit"])
        r_vec = e["r"] * vect_unit_abs
        # vitesse relative
        v_rel = v_translation + np.cross(omega, r_vec)
        V = np.linalg.norm(v_rel)
        if V < 1e-6:
            continue
        q = 0.5 * config.rho_air * V**2  # pression exercée par l'air sur le boomerang
        # calcul de la normale à l'élément de pale
        n = rot_current.apply(np.array([0, 0, 1]))
        dF_portance = q * e["dS"] * Cz_temp * n
        dF_trainee = q * e["dS"] * Cx_temp * (-v_rel / V)
        F_tot += dF_portance + dF_trainee
    return F_tot


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
