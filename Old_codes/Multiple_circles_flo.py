# Python program to Plot Circle

# importing libraries
import numpy as np
from matplotlib import pyplot as plt
import random
from scipy      import optimize
from matplotlib_scalebar.scalebar import ScaleBar
import pandas as pd
import scipy.stats as stats
from tkinter import filedialog
from tkinter.filedialog import asksaveasfilename
from sys import exit

# Creating equally spaced 100 data in range 0 to 2*pi
NumDrop = 100
x_fov = 20000
y_fov = 20000

#define primary variable
x_mo = np.array([None]*NumDrop)
y_mo = np.array([None]*NumDrop)
R_fit = np.array([None]*NumDrop)
center = np.array([None]*NumDrop)
Mol = np.array([None]*NumDrop)
Fl = np.array([None]*NumDrop)
centerf = np.zeros((NumDrop,2))
xc_2 = np.array([None]*NumDrop)
yc_2 = np.array([None]*NumDrop)
x_wr = np.array([None]*NumDrop)
y_wr = np.array([None]*NumDrop)
i = 0 
lower, upper = 0, 2*np.pi
# Generating x and y data
flag = 'false'
m_density = 0.1 ####Heterogeneous 0.05 and Homogeneous 0.1
while flag == 'false':
    
    radius = abs(np.random.normal(loc = 750  , scale =200))
    Nr = m_density*2*np.pi*radius
    sigma = radius/20
    mu_m, sigma_m = random.random()*2*np.pi, np.pi   ####Sigma_m  =  Heterogeneous (random.random() * (2.5 - 0.8) + 0.8) and Homogeneous np.pi
    Mol[i] = stats.truncnorm(
        (lower - mu_m) / sigma_m, (upper - mu_m) / sigma_m, loc=mu_m, scale=sigma_m)
    theta = Mol[i].rvs(int(Nr))
    x_mo[i] = radius * np.cos(theta) + np.random.randn(int(Nr))*sigma + random.randrange(0,x_fov)
    y_mo[i] = radius * np.sin(theta) + np.random.randn(int(Nr))*sigma + random.randrange(0,y_fov)
    xhe = np.hstack(x_mo[i])
    yhe = np.hstack(y_mo[i])
    F_density = 0.5
    radius_m = 3
    Nr_fl = F_density*2*np.pi*radius_m
    mu_f, sigma_f = random.random()*2*np.pi, np.pi  
    Fl[i] = stats.truncnorm(
        (lower - mu_f) / sigma_f, (upper - mu_f) / sigma_f, loc=mu_f, scale=sigma_f)
    sigma_fl = radius_m/0.2
    theta_fl = Fl[i].rvs(int(Nr_fl))
    x = np.array([None]*int(Nr))
    y = np.array([None]*int(Nr))
    
    kk = 0
    for kk in range(int(Nr)):
        x[kk] = radius_m * np.cos(theta_fl) + np.random.randn(int(Nr_fl))*sigma_fl + xhe[kk]
        y[kk] = radius_m * np.sin(theta_fl) + np.random.randn(int(Nr_fl))*sigma_fl + yhe[kk]
     
    xf = np.hstack(x) # final localisations
    yf = np.hstack(y)
    x_wr[i] = xf
    y_wr[i] = yf
    #Fitting to find the center in order to avoid cross section between droplets
    def calc_R(xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        return np.sqrt((xf-xc)**2 + (yf-yc)**2)

    def f_2(c):
        """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
        Ri = calc_R(*c)
        return Ri - Ri.mean()

    xm = np.mean(xf, axis = 0)
    ym = np.mean(yf, axis = 0)
    center_estimate = xm, ym
    center[i], ier = optimize.leastsq(f_2, center_estimate)
    centerf[i]= center[i].reshape(1,2)
    xc_2[i], yc_2[i] = centerf[i]
    Ri_fit       = calc_R(*centerf[i])
    R_fit[i]        = np.mean(Ri_fit)
    
    if i > 0:
        m = 1
        j = i
        while m==1 :
            R_cr = np.sqrt((centerf[i,0] - centerf[j - 1,0] )**2 + (centerf[i,1] - centerf[j - 1,1])**2)
            if R_cr < R_fit[i] + R_fit[j - 1]:
                i = i - 1
                m = 0;
            j = j - 1
            if j == 0:
              m = 0;          
    i = i + 1
    if i == NumDrop:
        flag = 'True'
       
#Data Extraction
for n in range (NumDrop):
    col1 = "x (nm)"
    col2 = "y (nm)"
    data = pd.DataFrame({col1:x_wr[n],col2:y_wr[n]})
    nn = n + 1
    #data.to_csv( r'C:\Users\JABER003\OneDrive - WageningenUR\Programming\Heterogeneity analyzer\Data_simulation\Droplet_%d.csv'%nn,  index=False)

col1 = "R (nm)"
data = pd.DataFrame({col1:R_fit})
data.to_csv( r'C:\Users\JABER003\OneDrive - WageningenUR\Programming\Heterogeneity analyzer\Data_simulation\R_fit.csv',  index=False)

# Plotting

plt.style.use('dark_background')
for i in range(NumDrop):
    plt.scatter(x_wr[i], y_wr[i] , s = 1, alpha = 0.1, color = 'c')
    ax=plt.gca()
    plt.xlim(0,x_fov)
    plt.ylim(0,y_fov)
    ax.set_aspect('equal')
    #plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
    plt.axis('off')
    #scalebar = ScaleBar(1, "nm", length_fraction=0.3, location = 'lower right', scale_loc = 'top',color = 'w', box_alpha = '0')
    ax.patch.set(alpha = 0.1)
    #ax.add_artist(scalebar)
plt.show()
a = asksaveasfilename(filetypes=(("PNG Image", "*.png"),("All Files", "*.*")), 
            defaultextension='.png')
if a:
    plt.savefig(a)  

#plt.savefig(r'C:\Users\JABER003\OneDrive - WageningenUR\Programming\Heterogeneity analyzer\Data_simulation\test.tif')
