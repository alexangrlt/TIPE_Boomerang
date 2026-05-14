import numpy as np
from scipy.interpolate import interp1d  # je vais en avoir besoin pr remplir le vide entre les valeurs que va me donner XFLR

def load_data(chemin_d_acces_au_fichier):
    data=np.loadtxt(chemin_d_acces_au_fichier,skiprows=............)
    
    #recup le alpha, Cl,Cd
    alpha=data[:,0]
    Cl=data[:,1]
    Cd=data[:,2]
    
    #puis utiliser interp1d dessus pour combler les trous
    Cl_p1d=interp1d(alpha,Cl,bounds_error=False,fill_value=(Cl[0],Cl[-1]))
    Cd_p1d=interp1d(alpha,Cd,bounds_error=False,fill_value=(Cd[0],Cd[-1]))
    
    return Cl_p1d,Cd_p1d    #ce sont des listes normalement