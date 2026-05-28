import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from core.boomerang_config import Boomerang_standard, BoomerangConfig
from core.blade_elements import get_blade_element, Cl_p1d, Cd_p1d


def simulate_projectile(position_init, vitesse_init, config, dt=0.0005, t_max=15):
    position = np.array(position_init, dtype=float)
    vitesse  = np.array(vitesse_init,  dtype=float)
    g = np.array([0.0, 0.0, -9.81])

    rot = R.from_euler('y', -80.0, degrees=True)

    I     = config.matrice_inertie()
    I_inv = np.linalg.inv(I)

    omega_monde = rot.apply(np.array([0.0, 0.0, 150.0]))

    elements = get_blade_element(config)

    pos_list    = []
    rot_list    = []
    omega_list  = []
    vit_list    = []   # DIAG : vitesse du CM
    F_aero_list = []   # DIAG : force aerodynamique totale
    incl_list   = []   # DIAG : inclinaison du plan (angle entre n_plan et vertical)
    t    = 0.0
    step = 0

    while t < t_max and position[2] >= 0.0:
        pos_list.append(position.copy())
        rot_list.append(rot.as_rotvec())
        omega_list.append(np.linalg.norm(omega_monde))
        vit_list.append(vitesse.copy())

        # Inclinaison : angle entre la normale au plan du boomerang et la verticale
        n_plan = rot.apply(np.array([0.0, 0.0, 1.0]))
        incl   = np.degrees(np.arccos(np.clip(abs(np.dot(n_plan, [0, 0, 1])), 0, 1)))
        incl_list.append(incl)

        F_aero, M_prec = compute_forces_be(elements, vitesse, omega_monde, rot, config)
        F_aero_list.append(F_aero.copy())
        F_tot = config.masse * g + F_aero

        rot_inv     = rot.inv()
        omega_corps = rot_inv.apply(omega_monde)
        M_corps     = rot_inv.apply(M_prec)

        def domega(oc, mc):
            return I_inv @ (mc - np.cross(oc, I @ oc))

        k1 = domega(omega_corps,             M_corps)
        k2 = domega(omega_corps + k1 * dt/2, M_corps)
        k3 = domega(omega_corps + k2 * dt/2, M_corps)
        k4 = domega(omega_corps + k3 * dt,   M_corps)
        omega_corps_new = omega_corps + (k1 + 2*k2 + 2*k3 + k4) * dt / 6

        TAU_SPIN = 5.0
        omega_corps_new[2] *= (1.0 - dt / TAU_SPIN)

        omega_monde_new = rot.apply(omega_corps_new)
        omega_moy = 0.5 * (omega_monde + omega_monde_new)
        rot = R.from_rotvec(omega_moy * dt) * rot

        omega_monde = omega_monde_new

        vitesse  += (F_tot / config.masse) * dt
        position += vitesse * dt

        step += 1
        # FIX bug 4 : renormalisation tous les 10 pas au lieu de 100
        # pour limiter la derive numerique du quaternion (omega~150 rad/s, dt=0.0005
        # => ~0.075 rad/pas, derive non negligeable sur 100 pas)
        if step % 10 == 0:
            q = rot.as_quat()
            rot = R.from_quat(q / np.linalg.norm(q))

        t += dt

    return {
        "px":     [p[0] for p in pos_list],
        "py":     [p[1] for p in pos_list],
        "pz":     [p[2] for p in pos_list],
        "pos":    pos_list,
        "rot":    np.array(rot_list),
        "omega":  np.array(omega_list),
        "vitesse": np.array(vit_list),
        "F_aero": np.array(F_aero_list),
        "incl":   np.array(incl_list),
        "dt":     dt,
    }


def compute_forces_be(elements, v_cm, omega_monde, rot, config):
    F_tot  = np.zeros(3)
    M_prec = np.zeros(3)

    n_plan = rot.apply(np.array([0.0, 0.0, 1.0]))

    for e in elements:
        axe_pale = rot.apply(e["vect_unit"])
        r_vec    = e["r"] * axe_pale

        v_rel = v_cm + np.cross(omega_monde, r_vec)
        V = np.linalg.norm(v_rel)
        if V < 0.5:
            continue

        v_rel_proj = v_rel - np.dot(v_rel, axe_pale) * axe_pale
        v_proj_mag = np.linalg.norm(v_rel_proj)
        if v_proj_mag < 1e-9:
            continue

        axe_corde_raw = np.cross(axe_pale, n_plan)
        norme_corde = np.linalg.norm(axe_corde_raw)
        if norme_corde < 0.1:
            axe_corde_raw = v_rel_proj - np.dot(v_rel_proj, n_plan) * n_plan
            norme_corde = np.linalg.norm(axe_corde_raw)
            if norme_corde < 1e-9:
                continue
        axe_corde = axe_corde_raw / norme_corde

        v_chordwise = np.dot(v_rel_proj, axe_corde)
        v_normal    = np.dot(v_rel_proj, n_plan)

        # FIX bug 3 : suppression du abs() sur v_chordwise pour conserver
        # le signe et distinguer la pale avancante de la pale en retraite.
        # Un v_chordwise negatif (pale en retraite) donne un alpha negatif
        # ou signe inverse, ce qui est physiquement correct.
        alpha_aero = np.degrees(np.arctan2(v_normal, v_chordwise + 1e-9))
        alpha = alpha_aero + np.degrees(e["twist"])

        Cl = float(Cl_p1d(alpha))
        Cd = float(Cd_p1d(alpha))
        q  = 0.5 * config.rho_air * V**2

        lift_dir = np.cross(v_rel_proj / v_proj_mag, axe_pale)
        ld_norm  = np.linalg.norm(lift_dir)
        if ld_norm < 1e-9:
            continue
        lift_dir = lift_dir / ld_norm

        # FIX bug 1 : la trainee est opposee a la vitesse relative DANS LE PLAN
        # de la section (v_rel_proj), pas a la vitesse 3D totale (v_rel).
        # Utiliser v_rel introduisait une composante axiale (le long de l'envergure)
        # inexistante en theorie BEM et qui perturbait le couple de precession.
        drag_dir = -v_rel_proj / v_proj_mag

        dF_lift = q * e["dS"] * Cl * lift_dir
        dF_drag = q * e["dS"] * Cd * drag_dir

        F_tot  += dF_lift + dF_drag
        # FIX bug 2 : le moment de precession inclut la trainee.
        # La trainee en bout de pale cree un couple gyroscopique non negligeable
        # qui contribue au virage. Ne compter que dF_lift sous-estimait la precession.
        M_prec += np.cross(r_vec, dF_lift + dF_drag)

    return F_tot, M_prec


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
