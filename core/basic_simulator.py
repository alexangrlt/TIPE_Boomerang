import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from core.boomerang_config import Boomerang_standard, BoomerangConfig
from core.blade_elements import get_blade_element, Cl_p1d, Cd_p1d


def simulate_projectile(position_init, vitesse_init, config, dt=0.0005, t_max=15):
    """
    Simule la trajectoire 3D d'un boomerang.

    Physique implementée :
    - Gravité
    - Forces aérodynamiques (BE method) avec déséquilibre pale avançante/reculante
    - Précession gyroscopique (eq. d'Euler dans le repere corps)
    - Rotation de l'orientation via quaternion

    Conventions :
    - Le boomerang est lancé dans la direction +X
    - Son plan est initialement quasi vertical (incliné ~80deg par rapport au sol)
    - omega est exprimé dans le REPERE MONDE
    - Les moments sont calculés dans le repere monde puis transformés en repere corps
      pour intégrer l'équation d'Euler
    """
    position = np.array(position_init, dtype=float)
    vitesse  = np.array(vitesse_init,  dtype=float)
    g = np.array([0.0, 0.0, -9.81])

    # Orientation initiale : plan du boomerang incliné a 80deg
    # Rotation autour de X de 80deg -> la normale au plan pointe vers le haut-avant
    rot = R.from_euler('x', 80.0, degrees=True)

    I     = config.matrice_inertie()       # tenseur inertie dans le repere corps
    I_inv = np.linalg.inv(I)

    # omega dans le repere monde, ~130 rad/s autour de la normale au plan
    omega_monde = rot.apply(np.array([0.0, 0.0, 130.0]))

    elements = get_blade_element(config)

    pos_list = []
    rot_list = []
    t = 0.0

    while t < t_max and position[2] >= 0.0:
        pos_list.append(position.copy())
        rot_list.append(rot.as_rotvec())

        # ---- Forces et moments dans le repere monde ----
        F_aero, M_monde = compute_forces_be(
            elements, vitesse, omega_monde, rot, config
        )
        F_tot = config.masse * g + F_aero

        # ---- Equation d'Euler dans le repere CORPS ----
        # On passe omega et M en repere corps pour l'integration
        rot_inv      = rot.inv()
        omega_corps  = rot_inv.apply(omega_monde)
        M_corps      = rot_inv.apply(M_monde)

        domega_corps = I_inv @ (M_corps - np.cross(omega_corps, I @ omega_corps))

        # Integration RK2 sur omega_corps
        omega_c2     = omega_corps + domega_corps * dt
        M_corps2     = M_corps  # approximation 1er ordre
        domega_c2    = I_inv @ (M_corps2 - np.cross(omega_c2, I @ omega_c2))
        omega_corps_new = omega_corps + 0.5 * (domega_corps + domega_c2) * dt

        # Clamp spin pour stabilite
        spd = np.linalg.norm(omega_corps_new)
        if spd > 250.0:
            omega_corps_new = omega_corps_new / spd * 250.0

        # Repasser omega en repere monde
        omega_monde = rot.apply(omega_corps_new)

        # ---- Integration translation ----
        vitesse  += (F_tot / config.masse) * dt
        position += vitesse * dt

        # ---- Mise a jour orientation ----
        # omega_monde * dt = vecteur de rotation infinitesimal
        rot = R.from_rotvec(omega_monde * dt) * rot

        t += dt

    return (
        [p[0] for p in pos_list],
        [p[1] for p in pos_list],
        [p[2] for p in pos_list],
        pos_list,
        np.array(rot_list),
    )


def compute_forces_be(elements, v_cm, omega_monde, rot, config):
    """
    Calcul des forces et moments aérodynamiques par la méthode des éléments de pale.

    Pour chaque tronçon :
    - Position r_vec = r * axe_pale (repere monde)
    - Vitesse locale = v_cm + omega x r_vec
    - alpha = angle entre v_rel et le plan du boomerang
    - Portance perpendiculaire a v_rel et a axe_pale
    - Trainee dans -v_rel
    """
    F_tot = np.zeros(3)
    M_tot = np.zeros(3)

    # Normale au plan du boomerang (repere monde)
    n_plan = rot.apply(np.array([0.0, 0.0, 1.0]))

    for e in elements:
        axe_pale = rot.apply(e["vect_unit"])  # repere monde
        r_vec    = e["r"] * axe_pale

        v_rel = v_cm + np.cross(omega_monde, r_vec)
        V = np.linalg.norm(v_rel)
        if V < 0.5:
            continue

        # Angle d'attaque : angle entre v_rel et le plan
        v_n     = np.dot(v_rel, n_plan)          # composante normale au plan
        v_t     = v_rel - v_n * n_plan           # composante dans le plan
        v_t_mag = np.linalg.norm(v_t)
        alpha   = np.degrees(np.arctan2(v_n, v_t_mag + 1e-9))

        Cl = float(Cl_p1d(alpha))
        Cd = float(Cd_p1d(alpha))

        q  = 0.5 * config.rho_air * V**2

        # Direction de portance : perp a v_rel ET perp a axe_pale
        # => c'est le vrai lift selon la theorie de l'aile
        lift_dir = np.cross(v_rel, axe_pale)
        ld_norm  = np.linalg.norm(lift_dir)
        if ld_norm < 1e-9:
            continue
        lift_dir = lift_dir / ld_norm

        # Signe : la portance s'oppose a la composante de v_rel selon n_plan
        if np.dot(lift_dir, n_plan) * v_n < 0:
            lift_dir = -lift_dir

        dF = q * e["dS"] * (Cl * lift_dir - Cd * v_rel / V)
        F_tot += dF
        M_tot += np.cross(r_vec, dF)

    return F_tot, M_tot


def plot_trajectory_3d(pos, title="Trajectoire du Boomerang"):
    pos = np.array(pos)
    fig = plt.figure()
    ax  = fig.add_subplot(projection="3d")
    ax.plot(pos[:, 0], pos[:, 1], pos[:, 2])
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")
    ax.set_zlabel("Z (m)")
    ax.set_title(title)
    plt.show()


def plot_rot(rotation, title="Orientation du Boomerang"):
    fig = plt.figure()
    ax  = fig.add_subplot(projection="3d")
    ax.plot(rotation[:, 0], rotation[:, 1], rotation[:, 2])
    ax.set_xlabel("X (rad)")
    ax.set_ylabel("Y (rad)")
    ax.set_zlabel("Z (rad)")
    ax.set_title(title)
    plt.show()


def plot_angles(rotation, dt):
    rotations = [R.from_rotvec(v) for v in rotation]
    angles    = np.array([r.as_euler("xyz", degrees=True) for r in rotations])
    t         = np.arange(len(rotation)) * dt
    plt.figure()
    plt.plot(t, angles[:, 0], label="Roulis (X)")
    plt.plot(t, angles[:, 1], label="Tangage (Y)")
    plt.plot(t, angles[:, 2], label="Lacet (Z)")
    plt.xlabel("Temps (s)")
    plt.ylabel("Angle (deg)")
    plt.legend()
    plt.grid(True)
    plt.show()
