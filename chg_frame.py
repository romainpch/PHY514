import numpy as np

def normalize(vector) :
    return vector / np.linalg.norm(vector)

def I2BF(Xeci, q) :
    #Passage de X dans le ECI à a dans le BF
    '''Body frame (cartesian) to Earth Centered Inertial frame (cartesian) thanks to the quaternion of the rotation'''
    x,y,z = Xeci
    a,b,c,d = q
    P = np.array([[2*(a**2+b**2)-1 , 2*(a*d+b*c) , 2*(b*d-a*c)],
                  [2*(b*c-a*d) , 2*(a**2+c**2)-1 , 2*(a*b+c*d)],
                  [2*(a*c+b*d) , 2*(c*d-a*b) , 2*(a**2+d**2)-1]])
    Xbf = np.dot(P,Xeci)
    return(Xbf)

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

