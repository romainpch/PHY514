import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import copy
import numpy as np
from datetime import date,datetime,timedelta,timezone
from skyfield import almanac
from astropy.io import fits
from math import *

from PIL import Image, ImageOps, ImageChops


fileradec = './globalstar_radec.csv'   
radec = pd.read_csv(fileradec)

radec['date'] = pd.to_datetime(radec['date'])
radec = radec[radec['len']>100] #On sélectionne les 4 strikes longs, les autres ne sont pas intéressants car redondants

#J'ai redéfini les valeurs de début et de fin de strike "à la main" car elles étaient pas bonnes dans le csv (étonnant d'ailleurs)
s_values = np.array([[2071.561768,1351.51206],
                     [1398.0,1462.6],
                     [1735.922266,1404.838922],
                     [1061.422021,1512.674988]])

e_values = np.array([[1735.8,1406.1],
                     [1060.891049,1512.340363],
                     [1400.152049,1458.731302],
                     [726.7,1567.0]])


def genere_light_curve(radec, N_strikes, s_values, e_values, i_plot = True, f_plot = True, pause_time=1., wide=5, inverse=True) :
    '''
    This function generates the light curve.
    It returns two arrays : the light curve itself and the array of times 

    radec       <class 'pandas.core.frame.DataFrame'>   List of all the strikes taken into account in the light curve generation
    N_strikes   <class 'int'>                           Number of strikes taken into account
    s_values    <class 'numpy.ndarray'>                 Array of the start coordinates of each strike, sorted as they are in radec
    e_values    <class 'numpy.ndarray'>                 Array of the end coordinates of each strike, sorted as they are in radec
    i_plot      <class 'bool'>                          If you want to inderstand what is going on with image rotation
    i_plot      <class 'bool'>                          If you want to plot the light curve
    pause_time  <class 'float'>                         Pause time used for each images
    wide        <class 'int'>                           Number of pixels summed
    inverse     <class 'bool'>                          Apply correction if the stike goes from right to left instead of left to right
    '''
    light_curve_global = np.array([])
    light_curve_times_global = np.array([])

    for i in range(N_strikes) :
        #selection of the corresponding .fits
        subradec = copy.deepcopy(radec[i:i+1])
        filefits = './18_34_35/'+os.path.basename(radec.iloc[i]['image'])

        hdu_list = fits.open(filefits)
        image_data = hdu_list[0].data
        hdu_list.close()

        dim_image = np.shape(image_data)

        
        S = np.array(s_values[i]).reshape((2,1))
        E = np.array(e_values[i]).reshape((2,1))
        s_x = S[0]
        s_y = S[1]
        e_x = E[0]
        e_y = E[1]

        #Rotation of the image : definition of the matrix and application of the good transformation
        angle_rad = atan((s_y-e_y)/(s_x-e_x))
        angle_deg = np.degrees(angle_rad)
        mat_rot = np.array([[cos(-angle_rad), -sin(-angle_rad)],
                            [sin(-angle_rad), cos(-angle_rad)]])

        shift = np.array([[int(dim_image[1]/2)],[int(dim_image[0]/2)]])

        S_rot = np.dot(mat_rot, S-shift) + shift
        E_rot = np.dot(mat_rot, E-shift) + shift

        image_data = Image.fromarray(image_data)
        image_data = np.array(image_data.rotate(angle_deg))

        tab = image_data[int(S_rot[1])-wide:int(S_rot[1])+wide , min(int(S_rot[0]),int(E_rot[0])):max(int(S_rot[0]),int(E_rot[0]))]
        light_curve = np.sum(tab,axis=0)
        light_curve_global = np.append(light_curve_global, light_curve)

        if i_plot :
            fig = plt.subplots(figsize=(15,8))


            plt.subplot(221)  # add subplot into first position in a 2x2 grid (upper left)
            plt.imshow(tab,cmap='gray')


            plt.subplot(223)  # add to third position in 2x2 grid (lower left)
            plt.plot([pause_time*i/len(light_curve) for i in range(len(light_curve))], light_curve, label='Light Curve')
            plt.ylabel('ADU')
            plt.xlabel('Temps (s)')


            plt.subplot(122) 
            plt.plot(S[0],S[1],'*',label='Start')
            plt.plot(E[0],E[1],'x',label='End')  
            plt.plot(S_rot[0],S_rot[1],'*',label='Start ROT')
            plt.plot(E_rot[0],E_rot[1],'x',label='End ROT')  
            plt.imshow(image_data)
            plt.legend()

        #Cette partie sert juste à se convaincre que les strikes sont espacés de 1 seconde et palier aux éventuels retards
        ref_time = np.array(['2021-01-25T18:34:42.784156000'], dtype='datetime64[ns]')
        s_time =  subradec['date']
        s_time = float((s_time.values - ref_time)*1.e-9) #la différence de temps est exprimée en nanosecondes d'où le *1.e-9 

        light_curve_times = np.linspace(s_time, s_time+pause_time, num=len(light_curve))
        light_curve_times_global = np.append(light_curve_times_global, light_curve_times)

    if inverse :
        light_curve_global = np.flip(light_curve_global)

    if f_plot :
        fig = plt.subplots(figsize=(15,8))
        plt.scatter(light_curve_times_global,light_curve_global,marker = '.', linewidths=0.3)
        plt.ylabel('ADU')
        plt.xlabel('Temps (s)')
    plt.show()

genere_light_curve(radec, 4, s_values, e_values, False)







#####################
# Poubelle à garder #
#####################

#s/e_values dans l'ordre naturel des strikes
# s_values = np.array([[1061.422021,1512.674988],
#                      [1398.0,1462.6],
#                      [1735.922266,1404.838922],
#                      [2071.561768,1351.51206]])

# e_values = np.array([[726.7,1567.0],
#                      [1060.891049,1512.340363],
#                      [1400.152049,1458.731302],
#                      [1735.8,1406.1]])