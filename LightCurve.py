from math import pi, log, sin, cos, exp
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import load, Topos, EarthSatellite
from astropy import units as u

def normalize(vector) :
    return vector / np.linalg.norm(vector)


####################################################
#Fonctions à coder / à verifier :
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

#Fonctions de conversion à vérifier, normalement c'est OK!
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
####################################################


class Surface :
    def __init__(self,vertices, normal, albedo=0.2) :
        self.P1 = vertices[0]
        self.P2 = vertices[1]
        self.P3 = vertices[2]
        self.normal_vector = normalize(normal)
        self.area = 0.5*np.linalg.norm(np.cross(self.P2-self.P1,self.P3-self.P1))
        self.albedo = albedo

class Object :
    def __init__(self,surf_list) :
        self.surf_list = surf_list

    def setObj(self,sat) :
        self.realobj = sat
    
    def setOrientation(self,orientation_quaternion) :
        self.orientation_quaternion = orientation_quaternion

    def computeBoxMesh(self, dim_x,dim_y,dim_z) :
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

    def computeSphereMesh(self, radius, N) :
        #a completer en prenant en compte N faces et non uniquement 264
        #il faut aussi changer les lignes l_surf.append(Surface en mettant des array + le vecteur normal pour aller avec la syntaxe de Surface
        l_surf = []
        top = 0,0,radius
        bottom = 0,0,-radius
        angle = np.pi/12
        
        m_vertex=[[(round(radius*np.sin(j*angle)*np.cos(2*i*angle),16),\
            round(radius*np.sin(j*angle)*np.sin(2*i*angle),16),\
            round(radius*np.cos(j*angle),16))\
            for i in range(13)] for j in range(1,12)]   #On cree tous les points de la triangulation

        for i in range(12): #Ajout des faces top et bottom
            l_surf.append(Surface(top,m_vertex[0][i+1],m_vertex[0][i]))
            l_surf.append(Surface(bottom,m_vertex[-1][i],m_vertex[-1][i+1]))
        
        for j in range(0,10):   #Ajout de toutes les faces latérales
            for i in range(12):
                l_surf.append(Surface(m_vertex[j][i],m_vertex[j][i+1],m_vertex[j+1][i]))
                l_surf.append(Surface(m_vertex[j][i+1],m_vertex[j+1][i+1],m_vertex[j+1][i]))
        self.surf_list = l_surf


class Sphere(Object) :
    def __init__(self,radius, nb_surf=264) :
        self.r = radius
        self.nb_surf = nb_surf
        Object.computeSphereMesh(self,r,N)


class Box(Object) :
    def __init__(self,dim_x,dim_y,dim_z) :
        self.nb_surf = 8
        Object.computeBoxMesh(self,dim_x,dim_y,dim_z)


def f(alpha, A_0=0.5, D=0.1, k=-0.5) :
    return A_0*exp(-alpha/D) + k*alpha + 1


#Les luminosités de Lommel-Seeliger et Lambert S_LS et S_L sont à coder
#Cependant elles fonctionnenet bien pour les planetes et astéoides
#Pour des satellites, préférer les modèles de Ashikhmin-Shirley et Cook-Torrance (a voir comment implémenter ça dans un second temps)
def uB(loc,t,q) :
    # pos_radec = loc.radec() #Obtenir ici les coordonées de loc sous forme np.array([r, ra, dec])
    pos_radec = loc.at(t).position.to(u.au).value
    return normalize(I2BF(radec2cart(pos_radec),q))

def S_LS(mu, mu_0) :
    # a coder, 
    return 1.0 

def S_L(mu, mu_0) :
    return 1.0

def S(mu,mu_0, alpha, c=0.1) :
    return f(alpha)*(S_LS(mu,mu_0) + c*S_L(mu,mu_0))

def L(object, obs, sun, t) :
    alpha = phase_angle(object.realobj, obs, t).to(u.deg).value
    res = 0
    for surf in object.surf_list :
        uB_n = surf.normal_vector
        uB_obs = uB(obs,t,object.orientation_quaternion)
        uB_sun = uB(sun,t,object.orientation_quaternion)

        mu = np.dot(uB_n, uB_obs)
        mu_0 = np.dot(uB_n, uB_sun)
        res += S(mu,mu_0, alpha)*surf.albedo*surf.area
    return res

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
# separations = phase_angle(iss,palaiseau,t)


#Génération d'une courbe de lumière :
satellite = Box(1.,1.,1.)
satellite.setObj(iss)
q_list = np.ones((len(t),4))

lightcurve = []
for i in range(len(t)) :
    q = q_list[i]
    satellite.setOrientation(q)
    lightcurve += [L(satellite,palaiseau,sun,t[i])]

plt.plot(t.tt,lightcurve)
plt.show()






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