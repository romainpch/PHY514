from math import pi, log, sin, cos, exp, acos
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, Topos, EarthSatellite
from astropy import units as u
from scipy.spatial.transform import Rotation as R
from scipy.signal import lombscargle
import random as rd

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
f_sample = 1000
alpha_deg = 95 #required to compare with Bradley and Axelrad results
satellite = BoxSatellite(3.,1.,1.)
#satellite = SphereSatellite(1.)
#satellite = CylinderSatellite(1.,1,1000)
w = 10*pi #satellite rotation

#Generating times and positions for the given phase angle
times = np.linspace(0,duration,num=duration*f_sample )
alpha_rad = pi*alpha_deg/180.
sat_pos = np.array([0.,0.,0.])
sun_pos =np.array([1.,0.,0.])
obs_pos =np.array([cos(alpha_rad),sin(alpha_rad),0])

alea = True #rot_axis is alea
Nb_genere = 1 #nb of generated lightcurves

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

# plt.xlabel('Rotational phase (deg)')
# plt.ylabel('Light curve')
# plt.show()



#########################
#  LightCurve analysis  #
#########################

def phase_folding(times,lightcurve,period,verbose=False): #Effectue le calcul des phases pour chaque point de mesure et les trie
    phases = [(t % period)/period for t in times]
    if verbose:
        print(phases)
    sorted_values = [val for _,val in sorted(zip(phases, lightcurve))]
    return sorted(phases),sorted_values

def phase_dispersion(sorted_values,mean=None): #Effectue le calcul de l'inverse du facteur de dispersion d'une courbe
    if mean is None:
        mean = np.mean(lightcurve)
    dispersion = 0
    var = 0
    n = len(sorted_values)
    sorted_values.append(sorted_values[0])
    for i in range(n):
        dispersion += (sorted_values[i]-sorted_values[i+1])**2
        var += (sorted_values[i]-mean)**2   #On pourrait pondérer en fonction de l'écart des phases
    return var/dispersion

def phase_reconstruction_diagram(times,lightcurve,periods): #Renvoie la liste des valeurs des facteurs de dispersion inverse correspondant à la liste de périodes en entrée
    mean = np.mean(lightcurve)
    phase_disp = [phase_dispersion(phase_folding(times,lightcurve,period)[1],mean) for period in periods]
    return phase_disp

def find_period(periods,phase_disp):    #Détermine la période pour laquelle l'inverse du facteur de dispersion est maximal
    t,m = max(zip(periods,phase_disp), key=(lambda x: x[1]))
    return t

def truncate(times,lightcurve,prop):    #Rééchantillone un jeu de données en choisissant aléatoirement une proportion donnée de ses valeurs
    n=len(times)
    sample = rd.sample(list(zip(times,lightcurve)),int(n*prop))
    sample.sort(key=(lambda x: x[0]))
    return [a for a,b in sample],[b for a,b in sample]


#Plot des courbes de lumière
plt.plot(times,lightcurve,color='r')
times_test,lightcurve_test = truncate(times,lightcurve,0.1)
lightcurve_test = [val+0.02*(rd.random()-1/2) for val in lightcurve_test]
plt.xlabel('Temps (s)')
plt.ylabel('Intensité lumineuse relative')
plt.plot(times_test,lightcurve_test,color='b')
plt.show()

#Calcul de la meileure période
log_periods = [k/100 for k in range(-200,0)]
periods = [10**l for l in log_periods]
freqs = [1/t for t in periods]
#Méthode de Lafler-Kinman
phase_disp = phase_reconstruction_diagram(times_test,lightcurve_test,periods)
prm, = plt.plot(freqs,normalize(phase_disp),color='b')
#Méthode de Lomb-Scargle
phase_disp_ls = lombscargle(times_test,lightcurve_test,[2*pi/t for t in periods])
ls, = plt.plot(freqs,normalize(phase_disp_ls),color='r')
#Combinaison
combine = [a*b for (a,b) in zip(phase_disp,phase_disp_ls)]
tot, = plt.plot(freqs,normalize(combine),color='g')
#Affichage
plt.xscale('log')
plt.xlabel('Fréquence (Hz)')
plt.legend([prm,ls,tot],['Lafler-Kinman','Lomb-Scargle','Combinaison'])
plt.show()

#Plot du diagramme de phase correspondant à la meilleure période
t = find_period(periods,combine)
plt.plot(*phase_folding(times_test,lightcurve_test,t),color='b')
plt.plot(*phase_folding(times,lightcurve,t),color='r')
plt.xlabel('Phase')
plt.ylabel('Intesité lumineuse relative')
plt.show()