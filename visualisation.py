import numpy as np
import math as m
from math import pi, cos, sin, acos
import matplotlib.pyplot as plt
#import skyfield


def dist(X) :
    '''Returns the EUclidian distance of X to (0,0,0)'''
    x,y,z = X
    return np.sqrt(x**2 + y**2 + z**2)

def sec2year(t) :
    return t/(3600*24*365)

def sec2days(t) :
    return t/(3600*24)

def cart2sph(X):
    x,y,z = X
    rho = m.sqrt(x**2 + y**2 + z**2)
    theta = m.acos(z/rho)
    phi = m.atan2(y,x)
    return((rho,theta%pi,phi%(2*pi)))

def sph2cart(X) :
    rho,theta,phi = X
    x = rho*sin(theta)*cos(phi)
    y = rho*sin(theta)*sin(phi)
    z = rho*cos(theta)
    return((x,y,z))

def sph2radec(X) :
    rho,theta,phi = X
    delta = pi/2-theta
    return((rho,phi,delta))

def Rot_quat(X, q) :
    '''Body frame (cartesian) to Earth Centered Inertial frame (cartesian) thanks to the quaternion of the rotation'''
    x,y,z = X
    X_vect = np.transpose(np.array([x,y,z]))
    a,b,c,d = q
    P= np.array([[2*(a**2+b**2)-1 , 2*(a*d+b*c) , 2*(b*d-a*c)],
                  [2*(b*c-a*d) , 2*(a**2+c**2)-1 , 2*(a*b+c*d)],
                  [2*(a*c+b*d) , 2*(c*d-a*b) , 2*(a**2+d**2)-1]])

    a = np.dot(X_vect,P)
    return((a[0],a[1],a[2]))

def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])

def cart2orbit(pos,vit):
    mu = 3.986*10**14
    energ = dist(vit)**2/2-mu/dist(pos)
    a = -mu/(2*energ)
    
    h = cross(pos,vit)
    i = m.acos(h[2]/dist(h))
    
    vh = cross(vit,h)
    evec = tuple(vh[j]/mu-pos[j]/dist(pos) for j in range(3))
    e = dist(evec)
    
    an = cross((0,0,1),h)
    raan = m.acos(an[0]/dist(an))
    raan = raan if an[1] > 0 else -raan
    
    arperi = m.acos((an[0]*evec[0]+an[1]*evec[1]+an[2]*evec[2])/(e*dist(an)))
    arperi = arperi if evec[2] > 0 else -arperi
    
    nu = m.acos((pos[0]*evec[0]+pos[1]*evec[1]+pos[2]*evec[2])/(e*dist(pos)))
    nu = nu if pos[0]*vit[0]+pos[1]*vit[1]+pos[2]*vit[2] > 0 else -nu
    
    return (a,i,e,raan,arperi,nu)


def print_ECO_ECI_radec(figname) :
    periods = []
    angular_velocity_ECI_ra = []
    angular_velocity_ECI_dec = []
    angular_velocity_ECO_lambda = []
    angular_velocity_ECO_theta = []

    for i in range(Nlines) :
        angular_velocity_BF = angular_velocities_BF[i]
        quaternion = orientations_quaternion[i]
        periods += [2*pi/dist(angular_velocity_BF)]

        angular_velocity_ECI_c = Rot_quat(angular_velocity_BF,quaternion)
        angular_velocity_ECI_radec = sph2radec(cart2sph(angular_velocity_ECI_c))
        angular_velocity_ECI_ra += [180.*angular_velocity_ECI_radec[1]/pi]
        angular_velocity_ECI_dec += [180.*angular_velocity_ECI_radec[2]/pi]

        _,inclination,_,Omega,_,_=cart2orbit(positions_ECI_c[i], speeds_ECI_c[i])

        angular_velocity_ECO_c = Rot_quat(Rot_quat(angular_velocity_ECI_c, (cos(Omega/2) , 0 , 0 , -sin(Omega/2))) , (cos(inclination/2) , -    sin(inclination/2) , 0 , 0 ) )
        angular_velocity_ECO_sph = cart2sph(angular_velocity_ECO_c)
        angular_velocity_ECO_theta += [180.*angular_velocity_ECO_sph[1]/pi] 
        angular_velocity_ECO_lambda += [180.*angular_velocity_ECO_sph[2]/pi]


    #Inertial frame figure
    fig1 = plt.figure(figsize=(10,3))

    plt.subplot(131)
    plt.plot(times, periods,linewidth=1.)
    plt.yscale('log')
    plt.xlabel('Time (year)')
    plt.ylabel('Period (s)')

    plt.subplot(132)
    plt.plot(times, angular_velocity_ECI_dec,linewidth=1.)
    plt.xlabel('Time (year)')
    plt.ylabel('Declination (deg)')
    # plt.xlim(0,8)
    plt.ylim(-90,10)

    plt.subplot(133)
    plt.plot(times, angular_velocity_ECI_ra,linewidth=1.)
    plt.xlabel('Time (year)')
    plt.ylabel('Right Ascension (deg)')
    # plt.xlim(0,8)
    plt.ylim(0,360)

    plt.tight_layout()
    plt.savefig(FigDestination + figname +'_inertial')

    #Orbital frame figure
    fig2 = plt.figure(figsize=(10,3))

    plt.subplot(131)
    plt.plot(times, angular_velocity_ECO_theta,linewidth=1.)
    plt.xlabel('Time (year)')
    plt.ylabel('$\\theta_{ECO}$ (deg)')
    # plt.xlim(0,8)

    plt.subplot(132)
    plt.plot(times, angular_velocity_ECO_lambda,linewidth=1.)
    plt.xlabel('Time (year)')
    plt.ylabel('$\lambda_{ECO}$ (deg)')
    # plt.xlim(0,8)
    # plt.ylim(0,360)

    ax = plt.subplot(133, projection='polar')
    ax.plot([ra*np.pi/180 for ra in angular_velocity_ECO_lambda],[180-dec for dec in angular_velocity_ECO_theta],linewidth=0.8)
    ax.set_rlabel_position(90)
    ax.set_rlim(0,90)
    ax.text(np.radians(90+10),ax.get_rmax()*3./4.,'$\\theta_{ECO}$',rotation=90,ha='center',va='center',size=6)
    ax.text(np.radians(20),ax.get_rmax()*1.1,'$\lambda_{ECO}$',rotation=-70,ha='center',va='center',size=6)
    ax.set_yticklabels(['165°','150°','135°','120°','105°'],size=6)

    plt.tight_layout()
    plt.savefig(FigDestination + figname +'_orbital')

def print_ECI_cart(figname,timescale,zoom=True) :
    fig1 = plt.figure(figsize=(10,8))

    plt.subplot(321)
    plt.plot(times, [180*a[0]/pi for a in angular_velocities_BF],linewidth=1.)
    plt.xlabel('Time ('+timescale+')')
    plt.ylabel('$\omega_x$')


    plt.subplot(323)
    plt.plot(times, [180*a[1]/pi for a in angular_velocities_BF],linewidth=1.)
    plt.xlabel('Time ('+timescale+')')
    plt.ylabel('$\omega_y$')


    plt.subplot(325)
    plt.plot(times, [180*a[2]/pi for a in angular_velocities_BF],linewidth=1.)
    plt.xlabel('Time ('+timescale+')')
    plt.ylabel('$\omega_z$')


    if zoom :
        plt.subplot(322)
        plt.plot(times[Nlines-80:], [180*a[0]/pi for a in angular_velocities_BF[Nlines-80:]],linewidth=1.)
        plt.xlabel('Time ('+timescale+')')
        plt.ylabel('$\omega_x$ (ZOOM)')

        plt.subplot(324)
        plt.plot(times[Nlines-80:], [180*a[1]/pi for a in angular_velocities_BF[Nlines-80:]],linewidth=1.)
        plt.xlabel('Time ('+timescale+')')
        plt.ylabel('$\omega_y$ (ZOOM)')

        plt.subplot(326)
        plt.plot(times[Nlines-870:], [ 180*a[2]/pi for a in angular_velocities_BF[Nlines-870:]],linewidth=1.)
        plt.xlabel('Time ('+timescale+')')
        plt.ylabel('$\Omega_z$ (ZOOM)')


    plt.tight_layout()
    plt.savefig(FigDestination + figname)





#Begining of the script
propa = 'propagation_v1.0.0_2020-10-22_11-02-54.txt' #LAGEOS-2 8-years propagation - step 1s
#propa = 'propagation_v1.0.0_2020-11-14_21-24-22.txt' #LAGEOS-2 8-years propagation - step 0.1s
#propa = 'propagation_v1.0.0_2020-12-09_22-31-07.txt' #Fig16 - OrtizGomez
#propa="propagation_v1.0.0_2020-12-11_09-51-13.txt" #Envisat -5 year

Path = './Simulations/LAGEOS-2/'
#Path = './Simulations/Envisat/Fig16_OrtizGomez/'
Nlines = 8*365
timescale='years'
FigDestination = './figures/LAGEOS-2/'


with open(Path + propa,'r') as data_file:
    lines = [line.strip('\n') for line in data_file.readlines()]

range_numbers = range(13,13+Nlines)

if timescale=='years' :
    times = [sec2year(float(lines[i].split()[0])) for i in range_numbers]
elif timescale=='days' :
    times = [sec2days(float(lines[i].split()[0])) for i in range_numbers]

speeds_ECI_c = [(float(lines[i].split()[1]),float(lines[i].split()[2]),float(lines[i].split()[3])) for i in range_numbers]
positions_ECI_c = [(float(lines[i].split()[4]),float(lines[i].split()[5]),float(lines[i].split()[6])) for i in range_numbers]
angular_velocities_BF = [(float(lines[i].split()[7]),float(lines[i].split()[8]),float(lines[i].split()[9])) for i in range_numbers]
orientations_quaternion = [(float(lines[i].split()[10]),float(lines[i].split()[11]),float(lines[i].split()[12]),float(lines[i].split()[13])) for i in range_numbers]

print_ECO_ECI_radec('lageos-2')
#print_ECI_cart('Envisat_Fig16_Velocity_cart.png','days')
#print_ECI_cart('Envisat_Velocity_cart.png',timescale,False)