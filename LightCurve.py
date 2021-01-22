from math import pi, log, sin, cos, exp, acos
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, Topos, EarthSatellite
from astropy import units as u


###########################################
#  Frame changements and basic functions  #
###########################################


def normalize(vector) :
    return vector / np.linalg.norm(vector)

def I2BF(X, q) :
    #Cette fonction qui vient de visualisation.py passe de X dans le BF à a dans le ECI, il faudrait l'inverser 
    '''Body frame (cartesian) to Earth Centered Inertial frame (cartesian) thanks to the quaternion of the rotation'''
    x,y,z = X
    X_vect = np.transpose(np.array([x,y,z]))
    a,b,c,d = q
    P= np.array([[2*(a**2+b**2)-1 , 2*(a*d+b*c) , 2*(b*d-a*c)],
                  [2*(b*c-a*d) , 2*(a**2+c**2)-1 , 2*(a*b+c*d)],
                  [2*(a*c+b*d) , 2*(c*d-a*b) , 2*(a**2+d**2)-1]])

    a = np.dot(X_vect,P)
    return((a[0],a[1],a[2]))

def sph2cart(Xsph) :
    rho, theta, phi = Xsph
    x = rho*sin(theta)*cos(phi)
    y = rho*sin(theta)*sin(phi)
    z = rho*cos(theta)
    return np.array([x,y,z])

def sph2radec(Xsph) :
    rho, theta, phi = Xsph
    delta = pi/2-theta
    return np.array([rho,phi,delta])

def radec2sph(Xradec) : 
    rho, phi, delta = Xradec
    theta = pi/2 - delta
    return np.array([rho, theta, phi])

def radec2cart(Xradec) :
    Xsph = radec2sph(Xradec)
    return sph2cart(Xsph)



##############################################
#  Definition of space objects and geometry  #
##############################################


class Surface :
    def __init__(self,vertices, normal_vector=None, albedo=0.2) :
        self.P1 = vertices[0]
        self.P2 = vertices[1]
        self.P3 = vertices[2]
        if normal_vector is None:
            self.normal_vector = normalize(np.cross(self.P2-self.P1,self.P3-self.P2))
        else : 
            self.normal_vector = normal_vector
        self.area = 0.5*np.linalg.norm(np.cross(self.P2-self.P1,self.P3-self.P1))
        self.albedo = albedo

class SpaceObject :
    def __init__(self,geometry=None,nb_surf=0) :
        self.geometry = geometry
        self.nb_surf = nb_surf
        self.position = np.array([0,0,0])

    def setObj(self,sat) :
        self.realobj = sat

    def computeSphereMesh(self,radius, N) :
        #a completer en prenant en compte N faces et non uniquement 264
        l_surf = []
        top = 0,0,radius
        bottom = 0,0,-radius
        angle = np.pi/12
        
        m_vertex=[[(round(radius*np.sin(j*angle)*np.cos(2*i*angle),16),\
            round(radius*np.sin(j*angle)*np.sin(2*i*angle),16),\
            round(radius*np.cos(j*angle),16))\
            for i in range(13)] for j in range(1,12)]   #On cree tous les points de la triangulation

        for i in range(12): #Ajout des faces top et bottom
            l_surf.append(Surface(np.array([top,m_vertex[0][i+1],m_vertex[0][i]])))
            l_surf.append(Surface(np.array([bottom,m_vertex[-1][i],m_vertex[-1][i+1]])))
        
        for j in range(0,10):   #Ajout de toutes les faces latérales
            for i in range(12):
                l_surf.append(Surface(np.array([m_vertex[j][i],m_vertex[j][i+1],m_vertex[j+1][i]])))
                l_surf.append(Surface(np.array([m_vertex[j][i+1],m_vertex[j+1][i+1],m_vertex[j+1][i]])))
        self.surf_list = l_surf

    def computeBoxMesh(self,dim_x,dim_y,dim_z) :
        x,y,z = dim_x/2, dim_y/2, dim_z/2
        l_surf = []
        #Top and bottom
        l_surf.append(Surface(np.array([[x,y,z],[-x,-y,z],[-x,y,z]]), np.array([0,0,1])))
        l_surf.append(Surface(np.array([[x,y,z],[x,-y,z],[-x,-y,z]]), np.array([0,0,1])))
        l_surf.append(Surface(np.array([[-x,-y,-z],[x,-y,-z],[-x,y,-z]]), np.array([0,0,-1])))
        l_surf.append(Surface(np.array([[x,-y,-z],[x,y,-z],[-x,y,-z]]), np.array([0,0,-1])))
        #Back and front
        l_surf.append(Surface(np.array([[x,y,z],[-x,y,z],[-x,y,-z]]), np.array([0,1,0])))
        l_surf.append(Surface(np.array([[x,y,z],[-x,y,-z],[x,y,-z]]), np.array([0,1,0])))
        l_surf.append(Surface(np.array([[-x,-y,z],[x,-y,z],[x,-y,-z]]), np.array([0,-1,0])))
        l_surf.append(Surface(np.array([[-x,-y,z],[x,-y,-z],[-x,-y,-z]]), np.array([0,-1,0])))
        #Left and right
        l_surf.append(Surface(np.array([[-x,y,z],[-x,-y,z],[-x,y,-z]]), np.array([-1,0,0])))
        l_surf.append(Surface(np.array([[-x,y,-z],[-x,-y,z],[-x,-y,-z]]), np.array([-1,0,0])))
        l_surf.append(Surface(np.array([[x,y,z],[x,y,-z],[x,-y,-z]]), np.array([1,0,0])))
        l_surf.append(Surface(np.array([[x,-y,z],[x,y,z],[x,-y,-z]]), np.array([1,0,0])))
        self.surf_list = l_surf  

class SphereSatellite(SpaceObject) :
    def __init__(self,radius, nb_surf=264) :
        SpaceObject.__init__(self,'sphere',nb_surf)
        SpaceObject.computeSphereMesh(self,radius,nb_surf)
        self.r = radius
        self.nb_surf = nb_surf
        
class BoxSatellite(SpaceObject) :
    def __init__(self,dim_x,dim_y,dim_z) :
        SpaceObject.__init__(self,'box', 12)
        SpaceObject.computeBoxMesh(self,dim_x,dim_y,dim_z)
        self.dimensions = [dim_x,dim_y,dim_z]

######################################
#  Functions to compute LightCurves  #
######################################


def angle_between(u,v) :
    return acos(np.dot(u,v))

def testreflectance(uB_sun,uB_n) :
    angle = angle_between(uB_sun,uB_n)
    if angle > pi/2 :
        return False
    if angle < 0 :
        return False
    return True

def phase_function(alpha, A_0=0.5, D=0.1, k=-0.5) :
    return A_0*exp(-alpha/D) + k*alpha + 1

#Les BRDF de Lommel-Seeliger et Lambert S_LS et S_L fonctionnent bien pour les planetes et astéoides
#Pour des satellites, préférer les modèles de Ashikhmin-Shirley et Cook-Torrance (a voir comment implémenter ça dans un second temps)
def brdf_function(mu,mu_0, alpha, c=0.1) :
    return phase_function(alpha)*mu*mu_0*(c + 1/(mu + mu_0))

def brightness(satellite, alpha,  observer_BF, sun_BF) :
    res = 0
    for surf in satellite.surf_list :
        uB_n = surf.normal_vector
        uB_obs = observer_BF #Because uB_obs is the vector from the observer to the target, and target is at [0,0,0] in the BF
        uB_sun = sun_BF
        
        isreflected = testreflectance(uB_sun,uB_n)
        
        if isreflected :
            mu = np.dot(uB_n, uB_obs)
            mu_0 = np.dot(uB_n, uB_sun)
            if mu*mu_0 != 0 :
                S = brdf_function(mu,mu_0, alpha)*surf.albedo*surf.area
                res += S
        
    return res

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
# Le satellite tourne selon l'axe z_BF (= z_IF) à une vitesse de 1 tour par seconde -> ω = 2π
# La durée d'observation est 10sec, la fréquence d'échantillonage est 100Hz

duration = 10
f_sample = 100
times = np.linspace(0,duration,num=duration*f_sample )

satellite = BoxSatellite(1.,1.,1.) #1m cubesat
omega = 2*pi
earth_phase = 5*pi/4
sun_phase = 3*pi/4
alpha = earth_phase-sun_phase

uBF_sun_ref = np.array([cos(sun_phase), sin(sun_phase), 0])

LightCurve = []
for t in times :
    uBF_sun = np.array([cos(omega*t + sun_phase), sin(omega*t + sun_phase), 0])
    uBF_earth = np.array([cos(omega*t + earth_phase), sin(omega*t + earth_phase), 0])

    LightCurve += [brightness(satellite,alpha,uBF_earth,uBF_sun)]


plt.plot(times[1:-1],LightCurve[1:-1],label='Light Curve')

plt.show()

# lightcurve = []
# for i in range(len(t)) :
#     q = q_list[i]
#     satellite.setOrientation(q)
#     geocentric_pos = satellite.realobj.at(t[i]).position.m
#     lightcurve += [magnitude(L(satellite,palaiseau,sun,t[i]),np.linalg.norm(geocentric_pos))]

# plt.plot(t.tt,lightcurve)
# plt.show()




# m_sun = -26.73 #Sun magnitude

# def p_diff(psi) :
#     return 2*(sin(psi) + (pi-psi)*cos(psi))/(3*pi)

# def m_obj(psi, d=1., R=6e5, rho_spec=0.175, rho_diff=0.175) :
#     return m_sun - 2.5*log(d**2 * (rho_spec/4 + rho_diff*p_diff(psi))/(R**2))

# def E_RSO(psi, d=1., R=6e5, rho_spec=0.175, rho_diff=0.175) :
#     return 5.6e10 * 10**(-0.4*m_obj(psi, d, R, rho_spec, rho_diff))

# psis = np.linspace(0,pi,100)
# ms = [m_obj(psi) for psi in psis]
# Es = [E_RSO(psi) for psi in psis]

# plt.subplot(211)
# plt.plot(psis,ms)
# plt.ylabel('magnitude')
# plt.subplot(212)
# plt.plot(psis,Es)
# plt.xlabel('angle (rad)')
# plt.ylabel('irradiance ($ph/s/m^2$)')

# plt.show()