"""
Configuration physique du boomerang 
"""


from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True) #TOUS LES SETTINGS SONT FROZEN/IMMUTABLES POUR EVITER LES ERREURS
class BoomerangConfig:
  """Settings Boomerang 2 pales (rectangulaires)"""
  
  #Parametres Geometriques
  L=0.4 #longueur de pale en m
  l=0.04 #largeur de pale en m
  h=0.004 #epaisseur de pale en m
  
  #Masse
  m=0.08 #masse en kg
  
  @property
  def envergure(self):
    """envergure du boomerang"""
    return self.L * np.sqrt(2) #sqrt(2) car on considère un carré et la diag d'un carré est la longueur *sqrt(2), en m


  def surface(self):
    """surface totale du boomerang"""
    return self.L * self.l *2 #calcul de surface, en m²