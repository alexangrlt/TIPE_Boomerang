import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from core.boomerang_config import Boomerang_standard, BoomerangConfig
from core.blade_elements import get_blade_element, Cl_p1d, Cd_p1d


def simulate_projectile(position_init, vitesse_init, config, dt=0.0005, t_max=15):
    position = np.array(position_init, dtype=float)
    vitesse  = np.array(vitesse_init,  dtype=float)
    g = np.array([0.0, 0.0, -9.81])

    # Inclinaison autour de Y : -80 deg (quasi-vertical, legerement plus dresse
    # qu'avant). En pratique un lanceur droitier incline le plan a ~80-85 deg
    # de l'horizontale pour favoriser la precession en debut de vol.
    rot = R.from_euler('y', -80.0, degrees=True)

    I     = config.matrice_inertie()
    I_inv = np.linalg.inv(I)

    # omega initial : ~1450 rpm = 150 rad/s autour de la normale au plan du boomerang
    omega_monde = rot.apply(np.array([0.0, 0.0, 150.0]))

    elements = get_blade_element(config)

    pos_list   = []
    rot_list   = []
    omega_list = []
    t    = 0.0
    step = 0

    while t < t_max and position[2] >= 0.0:
        pos_list.append(position.copy())
        rot_list.append(rot.as_rotvec())
        omega_list.append(np.linalg.norm(omega_monde))

        F_aero, M_prec = compute_forces_be(elements, vitesse, omega_monde, rot, config)
        F_tot = config.masse * g + F_aero

        # Equation d'Euler dans le repere corps
        rot_inv     = rot.inv()
        omega_corps = rot_inv.apply(omega_monde)
        M_corps     = rot_inv.apply(M_prec)

        def domega(oc, mc):
            return I_inv @ (mc - np.cross(oc, I @ oc))

        # RK4 sur omega_corps
        k1 = domega(omega_corps,             M_corps)
        k2 = domega(omega_corps + k1 * dt/2, M_corps)
        k3 = domega(omega_corps + k2 * dt/2, M_corps)
        k4 = domega(omega_corps + k3 * dt,   M_corps)
        omega_corps_new = omega_corps + (k1 + 2*k2 + 2*k3 + k4) * dt / 6

        # Freinage aerodynamique du spin.
        # TAU_SPIN = 5s : plus realiste pour un boomerang PLA leger (28g).
        # Un tau trop long (8s) maintient un spin fort mais ralentit la precession
        # en fin de vol car omega_prec = M / (I_zz * omega_spin) -> quand omega_spin
        # reste eleve, la precession reste lente si M_prec ne compense pas.
        # Avec tau=5s le spin chute plus vite, M_prec/I_zz/omega_spin augmente
        # en milieu de vol -> la boucle se ferme mieux.
        TAU_SPIN = 5.0
        omega_corps_new[2] *= (1.0 - dt / TAU_SPIN)

        # Integration rotation avec omega moyen (coherence avec RK4)
        omega_monde_new = rot.apply(omega_corps_new)
        omega_moy = 0.5 * (omega_monde + omega_monde_new)
        rot = R.from_rotvec(omega_moy * dt) * rot

        omega_monde = omega_monde_new

        # Integration translation
        vitesse  += (F_tot / config.masse) * dt
        position += vitesse * dt

        # Renormalisation quaternion toutes les 100 etapes
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
        np.array(omega_list),
    )


def compute_forces_be(elements, v_cm, omega_monde, rot, config):
    """
    Calcule les forces et moments aerodynamiques par la methode des elements de pale.

    CORRECTIONS APPLIQUEES :
    1. Angle d'attaque calcule dans le repere de section correct (via axe_corde).
    2. Le twist geometrique de chaque troncon est pris en compte dans alpha.
    3. Direction de portance corrigee : perp. a v_rel_proj dans le plan de section.
    4. Guard robuste sur axe_corde : quand le plan du boomerang devient quasi-
       horizontal (fin de vol), axe_corde peut devenir colineaire a axe_pale.
       On detecte ce cas et on skip le troncon proprement au lieu de diviser par zero.

    Retourne :
    - F_tot  : force totale (portance + trainee) pour la translation
    - M_prec : moment de precession (portance seule) pour la rotation
    """
    F_tot  = np.zeros(3)
    M_prec = np.zeros(3)

    # Normale au plan du boomerang dans le repere monde
    n_plan = rot.apply(np.array([0.0, 0.0, 1.0]))

    for e in elements:
        axe_pale = rot.apply(e["vect_unit"])
        r_vec    = e["r"] * axe_pale

        v_rel = v_cm + np.cross(omega_monde, r_vec)
        V = np.linalg.norm(v_rel)
        if V < 0.5:
            continue

        # Projection dans le plan de section (perp. a l'envergure)
        v_rel_proj = v_rel - np.dot(v_rel, axe_pale) * axe_pale
        v_proj_mag = np.linalg.norm(v_rel_proj)
        if v_proj_mag < 1e-9:
            continue

        # axe_corde : direction de corde locale = cross(axe_pale, n_plan)
        # Guard robuste : si axe_pale ~ n_plan (plan boomerang horizontal),
        # axe_corde devient quasi-nul. On utilise alors un axe de secours
        # base sur la vitesse relative projetee pour ne pas perdre la portance.
        axe_corde_raw = np.cross(axe_pale, n_plan)
        norme_corde = np.linalg.norm(axe_corde_raw)
        if norme_corde < 0.1:  # seuil plus robuste que 1e-9
            # Plan quasi-horizontal : la corde est dans le plan horizontal.
            # On prend la composante horizontale de v_rel_proj comme axe_corde.
            axe_corde_raw = v_rel_proj - np.dot(v_rel_proj, n_plan) * n_plan
            norme_corde = np.linalg.norm(axe_corde_raw)
            if norme_corde < 1e-9:
                continue
        axe_corde = axe_corde_raw / norme_corde

        # Composantes de v_rel_proj dans le repere de section
        v_chordwise = np.dot(v_rel_proj, axe_corde)
        v_normal    = np.dot(v_rel_proj, n_plan)

        # Angle d'attaque aerodynamique + twist geometrique local
        alpha_aero = np.degrees(np.arctan2(v_normal, abs(v_chordwise) + 1e-9))
        alpha = alpha_aero + np.degrees(e["twist"])  # twist negatif au bout -> reduit alpha

        Cl = float(Cl_p1d(alpha))
        Cd = float(Cd_p1d(alpha))
        q  = 0.5 * config.rho_air * V**2

        # Direction de portance : perp. a v_rel_proj ET a axe_pale
        lift_dir = np.cross(v_rel_proj / v_proj_mag, axe_pale)
        ld_norm  = np.linalg.norm(lift_dir)
        if ld_norm < 1e-9:
            continue
        lift_dir = lift_dir / ld_norm

        drag_dir = -v_rel / V

        dF_lift = q * e["dS"] * Cl * lift_dir
        dF_drag = q * e["dS"] * Cd * drag_dir

        F_tot  += dF_lift + dF_drag
        M_prec += np.cross(r_vec, dF_lift)

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
