from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


position = [0, 0, 10]  # x, y, z
vitesse = [5, 2, 0]

g = -9.81  # gravité en m/s^2
dt = 0.01  # pas de temps en secondes
t_tot = 5  # durée totale en secondes

Px = []
Py = []
Pz = []

t = 0

while t < t_tot and position[2] > 0:
    # Stockage AVANT calcul (pour avoir la position initiale)
    Px.append(position[0])
    Py.append(position[1])
    Pz.append(position[2])
    
    # Calculs physiques
    vitesse = [vitesse[0], vitesse[1], vitesse[2] + g*dt]
    position = [position[0] + vitesse[0]*dt, position[1] + vitesse[1]*dt, position[2] + vitesse[2]*dt]
    
    t += dt

# Plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot(Px, Py, Pz)
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
plt.show()
