# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 09:30:54 2020

@author: Maxence
"""

from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

file = 'D-SPOSE-master/input/sc_geometry.txt'

class Surface:  #Represente un triangle de la triangulation
    def __init__(self,_v1,_v2,_v3,_normal=None):
        if _normal is None:
            self.normal = give_normal(_v1,_v2,_v3)
        else:
            self.normal = _normal
        self.v1 = _v1
        self.v2 = _v2
        self.v3 = _v3

    def __repr__(self):
        normal_repr = repr(self.normal[0])+', '+repr(self.normal[1])+', '+repr(self.normal[2])
        v1_repr = repr(self.v1[0])+', '+repr(self.v1[1])+', '+repr(self.v1[2])
        v2_repr = repr(self.v2[0])+', '+repr(self.v2[1])+', '+repr(self.v2[2])
        v3_repr = repr(self.v3[0])+', '+repr(self.v3[1])+', '+repr(self.v3[2])
        return 'Normal : '+normal_repr +'\nVertex 1 : '+v1_repr+'\nVertex 2 : '+v2_repr+'\nVertex 3 :'+v3_repr+'\n'



def read_file(file):    #Pour lire un fichier de type 'sc_geometry.txt', renvoie une liste de surfaces
    l= []
    with open(file,'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0] != '%':
                values = line[:-1].split(sep='\t')
                normal = float(values[0]), float(values[1]), float(values[2])
                v1 = float(values[3]), float(values[4]), float(values[5])
                v2 = float(values[6]), float(values[7]), float(values[8])
                v3 = float(values[9]), float(values[10]), float(values[11])
                l.append(Surface(v1,v2,v3,normal))
    return(l)



def affichage(list_surf):   #Realise un affichage d'une triangulation
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    x,y,z = [],[],[]
    
    for surf in list_surf:
        x = [surf.v1[0], surf.v2[0],surf.v3[0]]
        y = [surf.v1[1], surf.v2[1],surf.v3[1]]
        z = [surf.v1[2], surf.v2[2],surf.v3[2]]

        #On effectue une rotation des points affichés, car il ne peut pas afficher deux points avec les même coordonnées X et Y
        c1,c2 = 0.999, 0.998
        s1,s2 = np.sqrt(1-c1**2),np.sqrt(1-c2**2)
        surface = np.array([[1,0,0],[0,c1,-s1],[0,s1,c1]]) @ np.array([[c2,0,s2],[0,1,0],[-s2,0,c2]]) @ np.array([x,y,z])
        
        ax.plot_trisurf(surface[0],surface[1],surface[2],linewidth=1,antialiased = True)
    plt.show()
    

            
def give_normal(v1,v2,v3):  #Renvoie le vecteur normal unitaire à un plan donne par trois points
    s1 = tuple(v2[i]-v1[i] for i in range(3))
    s2 = tuple(v3[i]-v1[i] for i in range(3))
    vect = s1[1]*s2[2]-s1[2]*s2[1],s1[2]*s2[0]-s1[0]*s2[2],s1[0]*s2[1]-s1[1]*s2[0]  #Produit scalaire
    return tuple(vect[i]/np.sqrt(vect[0]**2+vect[1]**2+vect[2]**2) for i in range(3))
      

          
def sphere(radius): #Renvoie la triangulation en 264 surfaces d'une sphere de rayon donne
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
    
    return l_surf



def box(dim_x,dim_y,dim_z): #Renvoie une triangulation d'un parallelepipede de dimensions donnees
    x,y,z = dim_x/2, dim_y/2, dim_z/2    
    l_surf = []
    
    #Top and bottom
    l_surf.append(Surface((x,y,z),(-x,-y,z),(-x,y,z)))
    l_surf.append(Surface((x,y,z),(x,-y,z),(-x,-y,z)))
    l_surf.append(Surface((-x,-y,-z),(x,-y,-z),(-x,y,-z)))
    l_surf.append(Surface((x,-y,-z),(x,y,-z),(-x,y,-z)))
    
    #Back and front
    l_surf.append(Surface((x,y,z),(-x,y,z),(-x,y,-z)))
    l_surf.append(Surface((x,y,z),(-x,y,-z),(x,y,-z)))
    l_surf.append(Surface((-x,-y,z),(x,-y,z),(x,-y,-z)))
    l_surf.append(Surface((-x,-y,z),(x,-y,-z),(-x,-y,-z)))
    
    #Left and right
    l_surf.append(Surface((-x,y,z),(-x,-y,z),(-x,y,-z)))
    l_surf.append(Surface((-x,y,-z),(-x,-y,z),(-x,-y,-z)))
    l_surf.append(Surface((x,y,z),(x,y,-z),(x,-y,-z)))
    l_surf.append(Surface((x,-y,z),(x,y,z),(x,-y,-z)))
    
    return l_surf



def export(list_surf,filename): #Ecrit un fichier pour une geometrie donnee au format D-SPOSE
    with open(filename,'w') as output:
        output.write('%%%%% Spacecraft Geometry\n')
        output.write('% Each row is one triangular surface\n')
        output.write('% Columns 1-3 are inward surface normal\n')
        output.write('% Columns 4-6, 7-9, 10-12 are coordinates of the vertices in body-fixed frame (m)\n')
        output.write('% Columns 13-18 are the coefficients of diffuse reflection, specular reflection, and absorption in (1) visible and (2) infrared spectrum\n')
        for surf in list_surf:
            for att in ['normal','v1','v2','v3']:
                for i in range(3):
                    output.write(str(format(getattr(surf,att)[i],'.16f')))
                    output.write('\t')
            output.write('0.3750000000000000\t0.2010000000000000\t0.4240000000000000\t0.095\t0.025\t0.880\n')   #Valeur arbitraire des coefficients de reflexion

affichage(box(10,5,26))