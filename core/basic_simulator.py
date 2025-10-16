import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R


def simulate_projectile(position_init, vitesse_init, dt=0.01, t_max=5):
  position = position_init.copy()
  vitesse = vitesse_init.copy()
  t = 0
  Px, Py, Pz = [], [], []
  R_rot,R_omega=[],[]
  g = -9.81
  #j'ai considéré le boomerang perpendiculaire par rapport au sol
  # donc lors du lancer il est tourné de 90deg selon x (avec x vers l'avant, y vers la gauche et z vers le haut) 
  # référentiel cartésien fixe dans le Référentiel terrestre (que je suppose à ce stade comme galiléen)
  rot_current=R.from_rotvec([np.pi/2,0,0])
  omega=np.array([0,8,0]) #on considère une vitesse angulaire constante ici
  
  while t < t_max and position[2] > 0:
    Px.append(position[0])
    Py.append(position[1])
    Pz.append(position[2])
    vitesse = [vitesse[0], vitesse[1], vitesse[2] + g*dt]
    position = [position[0] + vitesse[0]*dt, position[1] + vitesse[1]*dt, position[2] + vitesse[2]*dt]
    R_rot.append(rot_current.as_rotvec())
    rot_increment=R.from_rotvec(omega*dt)
    rot_current=rot_current*rot_increment
    t += dt
    
  rotation=np.array(R_rot)
  
  return Px, Py, Pz, rotation


def plot_rot(rotation, title="Position Angulaire"):
  fig = plt.figure()
  ax = fig.add_subplot(projection='3d')
  ax.plot(rotation[:,0], rotation[:,1],rotation[:,2], label='rotation')
  ax.set_xlabel('X (rad)')
  ax.set_ylabel('Y (rad)')
  ax.set_zlabel('Z (rad)')
  ax.set_title(title)
  ax.legend()
  plt.show()


def plot_trajectory_3d(Px, Py, Pz, title="Trajectoire du Projectile"):
  fig = plt.figure()
  ax = fig.add_subplot(projection='3d')
  ax.plot(Px, Py, Pz)
  ax.set_xlabel('X (m)')
  ax.set_ylabel('Y (m)')
  ax.set_zlabel('Z (m)')
  ax.set_title(title)
  plt.show()