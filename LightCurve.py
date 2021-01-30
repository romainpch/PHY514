from math import pi, log, sin, cos, exp, acos
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, Topos, EarthSatellite
from astropy import units as u
#custom libraries
from chg_frame import *
from SpaceObjectsGeometry import *

######################################
#  Functions to compute LightCurves  #
######################################


def phase_function(alpha, A_0=0.5, D=0.1, k=-0.5) :
    return A_0*exp(-alpha/D) + k*alpha + 1

#Les BRDF de Lommel-Seeliger et Lambert S_LS et S_L fonctionnent bien pour les planetes et astéoides
#Pour des satellites, préférer les modèles de Ashikhmin-Shirley et Cook-Torrance (a voir comment implémenter ça dans un second temps)
def brdf_function(mu,mu_0, alpha, c=0.1) :
    return phase_function(alpha)*mu*mu_0*(c + 1/(mu + mu_0))

def luminosity(sat,sat_pos,sun_pos,obs_pos,sat_att):
    sun_dir_BF = normalize(I2BF(sun_pos-sat_pos,sat_att))  #A vérifier que les changement de référentiel se font dans le bon sens
    obs_dir_BF = normalize(I2BF(obs_pos-sat_pos,sat_att))
    alpha = acos(np.dot(sun_dir_BF,obs_dir_BF))
    lum = 0
    for surf in sat.surf_list:
        surf_norm = surf.normal_vector
        mu = np.dot(surf_norm,obs_dir_BF)
        mu0 = np.dot(surf_norm,sun_dir_BF)
        if mu > 0 and mu0 > 0:
            lum += brdf_function(mu,mu0,alpha)*surf.albedo*surf.area
    return lum

def magnitude(L, d, m_ref = -26.7) :#référence par rapport au soleil
    '''Magnitude of an object of luminosity L at a distance d with respect to a reference object of magnitude m_ref'''
    E = L/(4*pi*d**2)
    return m_ref - 2.5*log(E)/log(10.)


####################
#  Skyfield part   #
####################


ts=load.timescale()
eph=load('de421.bsp')
sun, earth=eph['sun'],eph['earth']
palaiseau = Topos(latitude = 48.7107853, longitude = 2.2110581)

def phase_angle_pos(sat_pos,loc_pos,sol_pos):
    return (sol_pos-sat_pos).separation_from(loc_pos-sat_pos)

def phase_angle(sat,loc,t,sol=sun):
     sun_from_earth = sol-earth
     return phase_angle_pos(sat.at(t),loc.at(t),sun_from_earth.at(t))
    
iss = EarthSatellite('1 25544U 98067A   21011.39930122  .00001177  00000-0  29227-4 0  9996','2 25544  51.6469  38.2316 0000504 209.0007 284.1827 15.49278328264297',name='iss')
t=ts.utc(2021,1,11,16,range(30))


###########################
#  LightCurve generation  #
###########################


# On se donne un référentiel inertiel simplifié 
# Le satellite est un cube sat de 1m cube
# L'observateur, le soleil et le satellite sont fixes dans ce référentiel
# Ils ont pour positions Terre_IF(0,0,0)
#                        Soleil_IF(0,√2,0)
#                        Satellite_IF(√2/2,√2/2,0)
# L'angle de phase vaut toujours π/2
#
# Le BF du satelite est le référentiel inertiel translaté de telle sorte à ce que la position du satellite dans le BF soit (0,0,0)
# Ainsi, dans le BF, les positions sont à l'instant t : 
#                        Terre_BF(cos(ωt+5π/4),sin(ωt+5π/4),0)
#                        Soleil_BF(cos(ωt+3π/4),sin(ωt+3π/4),0)
#                        Satellite_BF(0,0,0)
#
# Le satellite tourne selon l'axe z_BF (= z_IF) à une vitesse de 1/4 tour par seconde -> ω = π/2
# La durée d'observation est 10sec, la fréquence d'échantillonage est 100Hz
            
duration = 10
f_sample = 100
times = np.linspace(0,duration,num=duration*f_sample )
satellite = BoxSatellite(1.,1.,1.)
sun_pos = np.array([0,0,0])
obs_pos = np.array([0,np.sqrt(2),0])
sat_pos = np.array([np.sqrt(2),np.sqrt(2),0])

rot_axis = normalize(np.array([0,0,1]))
w = pi/2
q_list = [(cos(w*t/2),rot_axis[0]*sin(w*t/2),rot_axis[1]*sin(w*t/2),rot_axis[2]*sin(w*t/2)) for t in times]

lightcurve = [luminosity(satellite,sat_pos,sun_pos,obs_pos,q) for q in q_list]
plt.plot(times,lightcurve,label='Light Curve')
plt.show()

