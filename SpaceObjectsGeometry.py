from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi

#Custom libraries
from chg_frame import normalize



class Surface :
    def __init__(self,vertices, normal_vector=None, albedo=0.2) :
        self.v1 = vertices[0]
        self.v2 = vertices[1]
        self.v3 = vertices[2]
        if normal_vector is None:
            self.normal_vector = normalize(np.cross(self.v2-self.v1,self.v3-self.v2))
        else : 
            self.normal_vector = normal_vector
        self.area = 0.5*np.linalg.norm(np.cross(self.v2-self.v1,self.v3-self.v1))
        self.albedo = albedo

    def __repr__(self):
            normal_vector_repr = repr(self.normal_vector[0])+', '+repr(self.normal_vector[1])+', '+repr(self.normal_vector[2])
            v1_repr = repr(self.v1[0])+', '+repr(self.v1[1])+', '+repr(self.v1[2])
            v2_repr = repr(self.v2[0])+', '+repr(self.v2[1])+', '+repr(self.v2[2])
            v3_repr = repr(self.v3[0])+', '+repr(self.v3[1])+', '+repr(self.v3[2])
            return 'Normal_vector : '+normal_vector_repr +'\nVertex 1 : '+v1_repr+'\nVertex 2 : '+v2_repr+'\nVertex 3 :'+v3_repr+'\n'


class SpaceObject :
    def __init__(self,geometry=None,nb_surf=0) :
        self.geometry = geometry
        self.nb_surf = nb_surf
        self.position = np.array([0,0,0])

    def setObj(self,sat) :
        self.realobj = sat

    def computeSphereMesh(self,radius, N_lat, N_long) :
        l_surf = []
        top = 0,0,radius
        bottom = 0,0,-radius
        angle_lat = np.pi/N_lat
        angle_long = np.pi/N_long
        
        m_vertex=[[(round(radius*np.sin(j*angle_lat)*np.cos(2*i*angle_long),16),\
            round(radius*np.sin(j*angle_lat)*np.sin(2*i*angle_long),16),\
            round(radius*np.cos(j*angle_lat),16))\
            for i in range(N_long+1)] for j in range(1,N_lat)]   #On cree tous les points de la triangulation

        for i in range(N_long): #Ajout des faces top et bottom
            l_surf.append(Surface(np.array([top,m_vertex[0][i+1],m_vertex[0][i]])))
            l_surf.append(Surface(np.array([bottom,m_vertex[-1][i],m_vertex[-1][i+1]])))
        
        for j in range(0,N_lat-2):   #Ajout de toutes les faces latérales
            for i in range(N_long):
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


    def computeCylinderMesh(self, radius, height, N):
        #Renvoie la triangulation en 4*N surfaces du cylindre de hauteur height et de rayon radius
        l_surf = []
        for k in range(N): 
            #top triangles :
            l_surf.append(Surface(np.array([[0.,0.,height],
                                            [radius*cos(2*k*pi/N), radius*sin(2*k*pi/N), height],
                                            [radius*cos(2*(k+1)*pi/N), radius*sin(2*(k+1)*pi/N), height]])))
            #bottom triangles :
            l_surf.append(Surface(np.array([[0.,0.,0.],
                                            [radius*cos(2*k*pi/N), radius*sin(2*k*pi/N), 0.],
                                            [radius*cos(2*(k+1)*pi/N), radius*sin(2*(k+1)*pi/N), 0.]])))    
            #face triangles
            l_surf.append(Surface(np.array([[radius*cos(2*k*pi/N), radius*sin(2*k*pi/N), 0.],
                                            [radius*cos(2*(k+1)*pi/N), radius*sin(2*(k+1)*pi/N), 0.],
                                            [radius*cos(2*k*pi/N),radius*sin(2*k*pi/N),height]])))

            l_surf.append(Surface(np.array([[radius*cos(2*k*pi/N),radius*sin(2*k*pi/N), height],
                                            [radius*cos(2*(k+1)*pi/N), radius*sin(2*(k+1)*pi/N), height],
                                            [radius*cos(2*(k+1)*pi/N), radius*sin(2*(k+1)*pi/N), 0.]]))) 
        self.surf_list = l_surf

    def affichage(self):   
        #Realise un affichage d'une triangulation
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        x,y,z = [],[],[]
        
        for surf in self.surf_list:
            x = [surf.v1[0], surf.v2[0],surf.v3[0]]
            y = [surf.v1[1], surf.v2[1],surf.v3[1]]
            z = [surf.v1[2], surf.v2[2],surf.v3[2]]

            #On effectue une rotation des points affichés, car il ne peut pas afficher deux points avec les même coordonnées X et Y
            c1,c2 = 0.999999999, 0.999999998
            s1,s2 = np.sqrt(1-c1**2),np.sqrt(1-c2**2)
            surface = np.array([[1,0,0],[0,c1,-s1],[0,s1,c1]]) @ np.array([[c2,0,s2],[0,1,0],[-s2,0,c2]]) @ np.array([x,y,z])
            
            
            ax.plot_trisurf(surface[0],surface[1],surface[2],linewidth=1,antialiased = True)
            #ax.auto_scale_xyz([-5, 5], [-5, 5], [-5, 5])
        plt.show()

    def export(self,filename): 
    #Ecrit un fichier pour une geometrie donnee au format D-SPOSE
        with open(filename,'w') as output:
            output.write('%%%%% Spacecraft Geometry\n')
            output.write('% Each row is one triangular surface\n')
            output.write('% Columns 1-3 are inward surface normal\n')
            output.write('% Columns 4-6, 7-9, 10-12 are coordinates of the vertices in body-fixed frame (m)\n')
            output.write('% Columns 13-18 are the coefficients of diffuse reflection, specular reflection, and absorption in (1) visible and (2) infrared spectrum\n')
            for surf in self.surf_list:
                for att in ['normal_vector','v1','v2','v3']:
                    for i in range(3):
                        output.write(str(format(getattr(surf,att)[i],'.16f')))
                        output.write('\t')
                output.write('0.3750000000000000\t0.2010000000000000\t0.4240000000000000\t0.095\t0.025\t0.880\n')   #Valeur arbitraire des coefficients de reflexion


############################
#  Specific Space Objects  #
############################

class SphereSatellite(SpaceObject) :
    def __init__(self,radius, n_lat=12, n_long=12) :
        SpaceObject.__init__(self,'sphere',(n_lat-1)*2*n_long)
        SpaceObject.computeSphereMesh(self,radius,n_lat,n_long)
        self.r = radius
        
class BoxSatellite(SpaceObject) :
    def __init__(self,dim_x,dim_y,dim_z) :
        SpaceObject.__init__(self,'box', 12)
        SpaceObject.computeBoxMesh(self,dim_x,dim_y,dim_z)
        self.dimensions = [dim_x,dim_y,dim_z]

class CylinderSatellite(SpaceObject) :
    def __init__(self,radius,height, nb_face_polygon=16) :
        SpaceObject.__init__(self,'box', 4*nb_face_polygon)
        SpaceObject.computeCylinderMesh(self,radius,height,nb_face_polygon)
        self.r = radius
        self.height = height

#####################
#  Other functions  #
#####################

def read_file(file):    
    #Pour lire un fichier de type 'sc_geometry.txt', renvoie une liste de surfaces
    surf_list = []
    with open(file,'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0] != '%':
                values = line[:-1].split(sep='\t')
                normal = float(values[0]), float(values[1]), float(values[2])
                v1 = float(values[3]), float(values[4]), float(values[5])
                v2 = float(values[6]), float(values[7]), float(values[8])
                v3 = float(values[9]), float(values[10]), float(values[11])
                surf_list.append(Surface(v1,v2,v3,normal))
    return(surf_list)
