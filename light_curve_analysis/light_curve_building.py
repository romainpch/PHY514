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

def plot_fits(filefits,radec,save=False):
    
  hdu_list = fits.open(filefits)
  image_data = hdu_list[0].data
  hdu_list.close()
  fig, ax1 = plt.subplots(figsize=(22,16))
  ax1.plot(radec['s_x'],radec['s_y'],'*',label='Start')
  ax1.plot(radec['e_x'],radec['e_y'],'x',label='End')    
  plt.imshow(image_data,cmap='gray')
  plt.colorbar()    
  plt.legend()
  if save==True:
    plt.savefig('ligh_curve')
  else:
    plt.show()

fileradec = './globalstar_radec.csv'   
radec = pd.read_csv(fileradec)
radec['date'] = pd.to_datetime(radec['date'])
radec = radec[radec['len']>50]
print(radec)
i = 3
subradec = copy.deepcopy(radec[i:i+1])
filefits = './18_34_35/'+os.path.basename(radec.iloc[i]['image'])
#plot_fits(filefits,subradec)


def lightcurve(filefits,radec) :
  hdu_list = fits.open(filefits)
  image_data = hdu_list[0].data
  hdu_list.close() 

  dim_image = np.shape(image_data)

  s_x = radec['s_x']
  s_y = radec['s_y']
  e_x = radec['e_x']
  e_y = radec['e_y']
  def  y(x) : #fonction affine passant par Start et End
    return s_y + (x - s_x)* (e_y-s_y)/(e_x-s_x)



  #Rotation 
  angle_rad = atan((s_y-e_y)/(s_x-e_x))
  angle_deg = np.degrees(angle_rad)
  mat_rot = np.array([[cos(-angle_rad), -sin(-angle_rad)],
                      [sin(-angle_rad), cos(-angle_rad)]])
  
  shift = np.array([[int(dim_image[1]/2)],[int(dim_image[0]/2)]])
  S = np.array([s_x,s_y])
  E = np.array([e_x,e_y])
  S_rot = np.dot(mat_rot, S-shift) + shift
  E_rot = np.dot(mat_rot, E-shift) + shift

  image_data = Image.fromarray(image_data)
  image_data = np.array(image_data.rotate(angle_deg))

  
  #fig = plt.subplots(figsize=(10,8))
  plt.subplot(211)
  wide = 5
  tab = image_data[int(S_rot[1])-wide:int(S_rot[1])+wide , min(int(S_rot[0]),int(E_rot[0])):max(int(S_rot[0]),int(E_rot[0]))]
  plt.imshow(tab,cmap='gray')
  plt.subplot(212)
  plt.plot(np.sum(tab,axis=0), label='Light Curve')
  plt.legend()


  fig = plt.subplots(figsize=(8,8))
  plt.plot(radec['s_x'],radec['s_y'],'*',label='Start')
  plt.plot(radec['e_x'],radec['e_y'],'x',label='End')  
  plt.plot(S_rot[0],S_rot[1],'*',label='Start ROT')
  plt.plot(E_rot[0],E_rot[1],'x',label='End ROT')  
  plt.imshow(image_data)
  plt.legend()
  plt.show()



fileradec = './globalstar_radec.csv'   
radec = pd.read_csv(fileradec)
radec['date'] = pd.to_datetime(radec['date'])
radec = radec[radec['len']>50]
print(radec)
i = 3
subradec = copy.deepcopy(radec[i:i+1])
filefits = './18_34_35/'+os.path.basename(radec.iloc[i]['image'])
lightcurve(filefits,subradec)






#Poubelle à garder : 

#generating light Curve
# nb_xpixels = abs(int(E_rot[0]-S_rot[0]))
# curve = np.zeros(nb_xpixels)
# coord_y = int(S_rot[1])
# for x in range(nb_xpixels) :
#   coord_x = x + int(S_rot[0])
#   lc_x = 0 
#   for i in range(2*wide +1) :
#     lc_x += image_data[coord_x, coord_y-wide+i]
#   curve[x] += lc_x
# plt.plot(np.flip(curve), label='generated curve')