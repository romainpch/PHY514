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

duration = 1
f_sample = 200
alpha_deg = 95 #required to compare with Bradley and Axelrad results
#satellite = BoxSatellite(1.,1.,1.)
#satellite = SphereSatellite(1.)
#satellite = CylinderSatellite(1.,1,30)
w = 2*pi #satellite rotation

#Generating times and positions for the given phase angle
times = np.linspace(0,duration,num=duration*f_sample )
alpha_rad = pi*alpha_deg/180.
sat_pos = np.array([0.,0.,0.])
sun_pos =np.array([1.,0.,0.])
obs_pos =np.array([cos(alpha_rad),sin(alpha_rad),0])

alea = False #rot_axis is alea
Nb_genere = 1 #nb of generated lightcurves

plt.title("Light curve for a phase angle : " + str(alpha_deg) + "deg")

if alea :
    for i in range(Nb_genere) :
        rot_axis = normalize(np.random.rand(3))
        q_list = [(cos(w*t/2),rot_axis[0]*sin(w*t/2),rot_axis[1]*sin(w*t/2),rot_axis[2]*sin(w*t/2)) for t in times]
        lightcurve = [luminosity(satellite,sat_pos,sun_pos,obs_pos,q) for q in q_list]
        plt.plot([360*t for t in times],lightcurve)
else :
    rot_axis = normalize(np.array([0,0,1]))
    q_list = [(cos(w*t/2),rot_axis[0]*sin(w*t/2),rot_axis[1]*sin(w*t/2),rot_axis[2]*sin(w*t/2)) for t in times]
    lightcurve = [luminosity(satellite,sat_pos,sun_pos,obs_pos,q) for q in q_list]
    plt.plot([360*t for t in times],lightcurve)

plt.xlabel('Rotational phase (deg)')
plt.ylabel('Light curve')
plt.show()