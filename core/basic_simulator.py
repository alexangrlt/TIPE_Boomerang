import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from core.boomerang_config import Boomerang_standard, BoomerangConfig
from core.blade_elements import get_blade_element, Cl_p1d, Cd_p1d


def simulate_projectile(position_init, vitesse_init, config, dt=0.0005, t_max=15):
    position = np.array(position_init, dtype=float)
    vitesse = np.array(vitesse_init, dtype=float)
    g = np.array([0.0, 0.0, -9.81])

    # Plan du boomerang incliné à 80deg par rapport au sol
    rot = R.from_euler("x", 80.0, degrees=True)

    I = config.matrice_inertie()
    I_inv = np.linalg.inv(I)

    # Rotation dans le sens anti-horaire vu du dessus (sens standard lancer droitier)
    # omega positif autour de la normale au plan = sens qui crée précession vers +Y
    omega_monde = rot.apply(np.array([0.0, 0.0, 130.0]))

    elements = get_blade_element(config)

    pos_list = []
    rot_list = []
    t = 0.0
    step = 0

    while t < t_max and position[2] >= 0.0:
        pos_list.append(position.copy())
        rot_list.append(rot.as_rotvec())

        F_aero, M_monde = compute_forces_be(elements, vitesse, omega_monde, rot, config)
        F_tot = config.masse * g + F_aero

        # Equation d'Euler dans le repère corps
        rot_inv = rot.inv()
        omega_corps = rot_inv.apply(omega_monde)
        M_corps = rot_inv.apply(M_monde)

        def domega(oc, mc):
            return I_inv @ (mc - np.cross(oc, I @ oc))

        # RK4 sur omega_corps
        k1 = domega(omega_corps, M_corps)
        k2 = domega(omega_corps + k1 * dt / 2, M_corps)
        k3 = domega(omega_corps + k2 * dt / 2, M_corps)
        k4 = domega(omega_corps + k3 * dt, M_corps)
        omega_corps_new = omega_corps + (k1 + 2 * k2 + 2 * k3 + k4) * dt / 6

        # Clamp spin
        spd = np.linalg.norm(omega_corps_new)
        if spd > 250.0:
            omega_corps_new = omega_corps_new / spd * 250.0

        omega_monde = rot.apply(omega_corps_new)

        # Integration translation
        vitesse += (F_tot / config.masse) * dt
        position += vitesse * dt

        # Mise à jour orientation
        rot = R.from_rotvec(omega_monde * dt) * rot

        # Renormalisation du quaternion toutes les 100 étapes
        step += 1
        if step % 100 == 0:
            q = rot.as_quat()
            rot = R.from_quat(q / np.linalg.norm(q))

        t += dt

    return (
        [p[0] for p in pos_list],
        [p[1] for p in pos_list],
        [p[2] for p in pos_list],
        pos_list,
        np.array(rot_list),
    )


def compute_forces_be(elements, v_cm, omega_monde, rot, config):
    F_tot = np.zeros(3)
    M_tot = np.zeros(3)

    # Normale au plan du boomerang (repère monde)
    # Pointe dans le sens du spin (convention droitier : vers haut-droite)
    n_plan = rot.apply(np.array([0.0, 0.0, 1.0]))

    for e in elements:
        axe_pale = rot.apply(e["vect_unit"])
        r_vec = e["r"] * axe_pale

        v_rel = v_cm + np.cross(omega_monde, r_vec)
        V = np.linalg.norm(v_rel)
        if V < 0.5:
            continue

        # Angle d'attaque
        v_n = np.dot(v_rel, n_plan)
        v_t = v_rel - v_n * n_plan
        v_t_mag = np.linalg.norm(v_t)
        alpha = np.degrees(np.arctan2(v_n, v_t_mag + 1e-9))

        Cl = float(Cl_p1d(alpha))
        Cd = float(Cd_p1d(alpha))
        q = 0.5 * config.rho_air * V**2

        # Portance : cross(v_rel, axe_pale), normalisé
        # Ce vecteur est perpendiculaire à v_rel ET à l'envergure de la pale
        lift_dir = np.cross(v_rel, axe_pale)
        ld_norm = np.linalg.norm(lift_dir)
        if ld_norm < 1e-9:
            continue
        lift_dir = lift_dir / ld_norm

        # La portance doit être du côté de n_plan (face aspirée du profil)
        if np.dot(lift_dir, n_plan) < 0:
            lift_dir = -lift_dir

        dF = q * e["dS"] * (Cl * lift_dir - Cd * v_rel / V)
        F_tot += dF
        M_tot += np.cross(r_vec, dF)

    return F_tot, M_tot


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


def plot_rot(rotation, title="Orientation du Boomerang"):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.plot(rotation[:, 0], rotation[:, 1], rotation[:, 2])
    ax.set_xlabel("X (rad)")
    ax.set_ylabel("Y (rad)")
    ax.set_zlabel("Z (rad)")
    ax.set_title(title)
    plt.show()


def plot_angles(rotation, dt):
    rotations = [R.from_rotvec(v) for v in rotation]
    angles = np.array([r.as_euler("xyz", degrees=True) for r in rotations])
    t = np.arange(len(rotation)) * dt
    plt.figure()
    plt.plot(t, angles[:, 0], label="Roulis (X)")
    plt.plot(t, angles[:, 1], label="Tangage (Y)")
    plt.plot(t, angles[:, 2], label="Lacet (Z)")
    plt.xlabel("Temps (s)")
    plt.ylabel("Angle (deg)")
    plt.legend()
    plt.grid(True)
    plt.show()
