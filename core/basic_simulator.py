import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import (
    interp1d,
)  # je vais en avoir besoin pr remplir le vide entre les valeurs que va me donner XFLR
from core.boomerang_config import Boomerang_standard, BoomerangConfig
from core.blade_elements import get_blade_element, Cl_p1d, Cd_p1d
from core.xflr_data import load_data

#! equation d'euler pour l'évolution de omega pr l'effet gyroscopique


def simulate_projectile(position_init, vitesse_init, config, dt=0.0001, t_max=20):
    position = position_init.copy()
    vitesse = vitesse_init.copy()
    t = 0
    Px, Py, Pz = [], [], []
    pos = []
    R_rot, R_omega = [], []
    g = 9.81

    # ? j'ai considéré le boomerang a 25° p/r a la normale au sol (z)
    # ? donc 65° p/r a x (avec x vers l'avant, y vers la gauche et z vers le haut)

    rot_current = R.from_rotvec([20 * np.pi / 180, 0, 0])
    # omega = np.array([0, 8, 0])  # on considère une vitesse angulaire constante ici

    elements = get_blade_element(config)

    I=config.matrice_inertie()
    I_inv=np.linalg.inv(I)
    omega=np.array([0.0, 0.0, 50.0])    #omega est mtn une variable contrairement a ce que j'avais fait avant avec la ligne 26

    while t < t_max and position[2] > 0:
        pos.append(
            [position[0], position[1], position[2]]
        )  # liste de listes des positions
        Px.append(position[0])
        Py.append(position[1])
        Pz.append(position[2])
        """Forces"""
        F_gravite = np.array([0, 0, -config.masse * g])

        F_pale, M_aero = compute_forces_be(elements, vitesse, omega, rot_current, config)

        F_tot = F_gravite + F_pale
        
        #? mise en place precession gyrosc
        gyro=np.cross(omega,I@omega)
        omega_point=I_inv@(M_aero-gyro)
        omega+=omega_point*dt

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

    return Px, Py, Pz, pos, rotation


def compute_forces_be(elements, v_translation, omega, rot_current, config):
    """calcule la force totale F_tot exercée sur chaque éléments

    Args:
        elements (dict): dictionnaire ayant toutes les infos de chaque element du boomerang
        v_translation (array): tableau des vitesses de translation du centre de masse du boomerang en fonction du temps en 3axes
        omega (array): tableau des vitesses de rotation 3axes en fonction du temps
        rot_current (scipy rot): orientation du boomerang
        config (_type_): appel des données de config
    """
    F_tot = np.zeros(3)
    M_tot = np.zeros(3)

    for e in elements:
        # position du troncon dans le repère terrestre (absolu)
        vect_unit_abs = rot_current.apply(e["vect_unit"])
        r_vec = e["r"] * vect_unit_abs
        # vitesse relative
        v_rel = v_translation + np.cross(omega, r_vec)
        V = np.linalg.norm(v_rel)
        
        if V < 1e-6:
            continue
        
        q = 0.5 * config.rho_air * V**2  # pression dynamique

        # calcul de la normale à l'élément de pale
        n = rot_current.apply(np.array([0, 0, 1]))
        if n[2] < 0:
            n = -n

        # angle d'attaque local : angle entre v_rel et le plan de pale
        v_normale = np.dot(v_rel, n)
        v_tang = np.linalg.norm(v_rel - v_normale * n)
        alpha_local = np.degrees(np.arctan2(v_normale, v_tang + 1e-9))

        # Cl et Cd depuis les polaires XFLR5
        Cl_temp = float(Cl_p1d(alpha_local))   # portance
        Cd_temp = float(Cd_p1d(alpha_local))   # traînée

        dF_portance = q * e["dS"] * Cl_temp * n
        dF_trainee  = q * e["dS"] * Cd_temp * (-v_rel / V)
        dF = dF_portance + dF_trainee
        F_tot += dF
        dM = np.cross(r_vec, dF)
        M_tot += dM

    return F_tot, M_tot


def plot_rot(rotation, title="Position Angulaire du centre d'inertie du Boomerang"):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot(rotation[:, 0], rotation[:, 1], rotation[:, 2], label="rotation")
    ax.set_xlabel("X (rad)")
    ax.set_ylabel("Y (rad)")
    ax.set_zlabel("Z (rad)")
    ax.set_title(title)
    ax.legend()
    plt.show()


def plot_trajectory_3d(pos, title="Trajectoire du Boomerang"):
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
    rotations = [R.from_rotvec(vec) for vec in rotation]
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
