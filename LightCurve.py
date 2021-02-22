from math import pi, log, sin, cos, exp, acos
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, Topos, EarthSatellite
from astropy import units as u
from scipy.spatial.transform import Rotation as R
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
    sun_dir_BF = normalize(change_ref(sun_pos-sat_pos,sat_att))  #A vérifier que les changement de référentiel se font dans le bon sens
    obs_dir_BF = normalize(change_ref(obs_pos-sat_pos,sat_att))
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


##########################
# Generating orientation #
##########################


def generate_attitude_list(rot_init, rot_vect, n_iter):
    rot_list = [rot_init]*n_iter
    rot_cur = rot_init
    for i in range(n_iter):
        rot_list[i] = rot_cur
        rot_cur = rot_vect * rot_cur
    return rot_list


###########################################
#  LightCurve generation : Basic example  #
###########################################

duration = 1
f_sample = 200
alpha_deg = 95 #required to compare with Bradley and Axelrad results
"""
satellite = BoxSatellite(3.,1.,1.)
#satellite = SphereSatellite(1.)
#satellite = CylinderSatellite(1.,1,1000)
w = 2*pi #satellite rotation

#Generating times and positions for the given phase angle
times = np.linspace(0,duration,num=duration*f_sample )
alpha_rad = pi*alpha_deg/180.
sat_pos = np.array([0.,0.,0.])
sun_pos =np.array([1.,0.,0.])
obs_pos =np.array([cos(alpha_rad),sin(alpha_rad),0])

alea = False #rot_axis is alea
Nb_genere = 20 #nb of generated lightcurves

#plt.title("Light curve for Sphere Saatellite at phase angle : " + str(alpha_deg) + "deg")

if alea :
    for i in range(Nb_genere) :
        rot_axis = normalize(np.random.rand(3))
        q_list = generate_attitude_list(R.identity(),R.from_rotvec(w/f_sample*rot_axis),duration*f_sample)
        lightcurve = [luminosity(satellite,sat_pos,sun_pos,obs_pos,q) for q in q_list]
        plt.plot([360*t for t in times],lightcurve)
else :
    rot_axis = normalize(np.array([0,0,1]))
    q_list = generate_attitude_list(R.identity(),R.from_rotvec(w/f_sample*rot_axis),duration*f_sample)
    lightcurve = [luminosity(satellite,sat_pos,sun_pos,obs_pos,q) for q in q_list]
    plt.plot([360*t for t in times],lightcurve)

plt.xlabel('Rotational phase (deg)')
plt.ylabel('Light curve')
plt.show()
"""


#########################
#  LightCurve analysis  #
#########################

def phase_folding(times,lightcurve,period,verbose=False):
    phases = [(t % period)/period for t in times]
    if verbose:
        print(phases)
    sorted_values = [val for _,val in sorted(zip(phases, lightcurve))]
    return sorted(phases),sorted_values

def phase_dispersion(sorted_values,mean=None):
    if mean is None:
        mean = np.mean(lightcurve)
    dispersion = 0
    var = 0
    n = len(sorted_values)
    sorted_values.append(sorted_values[0])
    for i in range(n):
        dispersion += (sorted_values[i]-sorted_values[i+1])**2
        var += (sorted_values[i]-mean)**2
    return var/dispersion
        
def phase_reconstruction_diagram(times,lightcurve,periods):
    mean = np.mean(lightcurve)
    phase_disp = [phase_dispersion(phase_folding(times,lightcurve,period)[1],mean) for period in periods]
    return phase_disp

def find_period(periods,phase_disp):
    t,m = max(zip(periods,phase_disp), key=(lambda x: x[1]))
    return t

duration = 5
f_sample = 200
times = np.linspace(0,duration,num =duration*f_sample)
print(len(times))
w = 10*pi
lightcurve = [np.sin(w*t)**3+0.4*np.random.random()+0.3*np.sin(0.2*pi*t) for t in times]
plt.plot(times,lightcurve)
plt.show()

log_periods = [k/1000 for k in range(-1000,0)]
periods = [10**l for l in log_periods]
phase_disp = phase_reconstruction_diagram(times,lightcurve,periods)
plt.plot(periods,phase_disp)
plt.xscale('log')
plt.show()

t = find_period(periods,phase_disp)
print(t)
plt.plot(*phase_folding(times,lightcurve,t,verbose=False))
plt.show()

    
    


