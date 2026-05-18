import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from scipy.interpolate import (
    interp1d,
)
from core.boomerang_config import Boomerang_standard, BoomerangConfig
from core.blade_elements import get_blade_element, Cl_p1d, Cd_p1d
from core.xflr_data import load_data

#! equation d'euler pour l'évolution de omega pr l'effet gyroscopique


def simulate_projectile(position_init, vitesse_init, config, dt=0.0005, t_max=15):
    position = np.array(position_init, dtype=float)
    vitesse = np.array(vitesse_init, dtype=float)
    t = 0
    Px, Py, Pz = [], [], []
    pos = []
    R_rot = []
    g = 9.81

    # Boomerang lancé quasi vertical, incliné ~20° sur X (axe de lancer)
    # Le plan du boomerang est presque perpendiculaire au sol au moment du lancer
    rot_current = R.from_rotvec([20 * np.pi / 180, 0, 0])

    elements = get_blade_element(config)

    I = config.matrice_inertie()
    I_inv = np.linalg.inv(I)

    # omega dans le repère monde : ~70 rad/s ≈ 11 tours/s, axe Z (rotation principale)
    # On tourne d'abord dans le repère boomerang puis on exprime dans le monde
    omega = rot_current.apply(np.array([0.0, 0.0, 70.0]))

    while t < t_max and position[2] > 0:
        pos.append([position[0], position[1], position[2]])
        Px.append(position[0])
        Py.append(position[1])
        Pz.append(position[2])

        F_gravite = np.array([0, 0, -config.masse * g])

        F_pale, M_aero = compute_forces_be(elements, vitesse, omega, rot_current, config)

        F_tot = F_gravite + F_pale

        # Précession gyroscopique (équation d'Euler)
        gyro = np.cross(omega, I @ omega)
        
        def omega_deriv(w):
            return I_inv @ (M_aero - np.cross(w, I @ w))

        k1 = omega_deriv(omega)
        k2 = omega_deriv(omega + k1 * dt)
        omega = omega + 0.5 * (k1 + k2) * dt

        # Limiter omega pour éviter divergence numérique
        omega_norm = np.linalg.norm(omega)
        if omega_norm > 200.0:
            omega = omega / omega_norm * 200.0

        acceleration = F_tot / config.masse
        vitesse += acceleration * dt

        position = position + vitesse * dt

        R_rot.append(rot_current.as_rotvec())
        rot_increment = R.from_rotvec(omega * dt)
        rot_current = rot_current * rot_increment
        t += dt

    rotation = np.array(R_rot)
    return Px, Py, Pz, pos, rotation


def compute_forces_be(elements, v_translation, omega, rot_current, config):
    """Calcule la force totale et le moment aérodynamique sur toutes les pales."""
    F_tot = np.zeros(3)
    M_tot = np.zeros(3)

    for e in elements:
        vect_unit_abs = rot_current.apply(e["vect_unit"])
        r_vec = e["r"] * vect_unit_abs

        v_rel = v_translation + np.cross(omega, r_vec)
        V = np.linalg.norm(v_rel)

        if V < 1e-6:
            continue

        q = 0.5 * config.rho_air * V**2

        # Normale au plan de pale (axe Z du boomerang dans le repère monde)
        n = rot_current.apply(np.array([0, 0, 1]))
        if n[2] < 0:
            n = -n

        # Angle d'attaque local
        v_normale = np.dot(v_rel, n)
        v_tang = np.linalg.norm(v_rel - v_normale * n)
        alpha_local = np.degrees(np.arctan2(v_normale, v_tang + 1e-9))

        Cl_temp = float(Cl_p1d(alpha_local))
        Cd_temp = float(Cd_p1d(alpha_local))

        dF_portance = q * e["dS"] * Cl_temp * n
        dF_trainee  = q * e["dS"] * Cd_temp * (-v_rel / V)
        dF = dF_portance + dF_trainee
        F_tot += dF
        M_tot += np.cross(r_vec, dF)

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
