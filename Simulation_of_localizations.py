# importing libraries
import sys
import os
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
import random
import pandas as pd
from tkinter import filedialog
import tkinter as tk
#Protein density
m_density = 0.6
#Number of droplets
NumDrop = 100
#Dimension of the field of view in nanometer
x_fov = 30000
y_fov = 30000
#Size of the protein in nm
d_pr = 22
#Localization accuracy in nanometer
sigma = 30
#Aggregation probablity term
agg = 0.0
#Mean valu for the exponential distribution of the number of localizations per molecule
meann = 5
#Localization Precision in nanometer
sigma_fl = 50

x_mo = np.array([None]*NumDrop)
R_simu = np.array([None]*NumDrop)
y_mo = np.array([None]*NumDrop) 
x_mo2 = np.array([None]*NumDrop)
y_mo2 = np.array([None]*NumDrop)
x_wr = np.array([None]*NumDrop)
y_wr = np.array([None]*NumDrop)
center = np.array([None]*NumDrop)
centerf = np.zeros((NumDrop,2))
xc_2 = np.array([None]*NumDrop)
yc_2 = np.array([None]*NumDrop)
R_fit = np.array([None]*NumDrop)
#Counter for the number of proteins
co_mul = 0
#for droplet counting
flag = 'false'
while flag == 'false':
    # Setting for droplets' radius
    radius = abs(np.random.normal(loc = 750  , scale =200))
    #Counter for independant proteins in aggregation mode (agg is not equal to 0) 
    ind = 0
    #Number of proteins
    Nr = int(m_density*2*np.pi*radius/d_pr)
    randr1 = np.array([None]*Nr)
    randr2 = np.array([None]*Nr)

    randr1 = np.random.normal(loc = sigma, scale = 0, size = Nr)
    randr2 = np.random.normal(loc = sigma, scale = 0, size = Nr)        
    theta = np.random.permutation(np.linspace(0,2*np.pi,Nr)) 
    light = 'red'
    cts =1
    x_cts = np.array([None]*Nr)
    L = np.array([None]*Nr)
    y_cts = np.array([None]*Nr)
    agg_co = np.zeros(Nr)
    ran = np.zeros(int(Nr))
    ran_p = np.zeros(int(Nr))
    ran_n = np.zeros(int(Nr))
    # Generating x and y molecule position
    x_A = random.randrange(0,x_fov)
    y_A = random.randrange(0,y_fov)
    x_mo[co_mul] = (radius + (randr1)) * np.cos(theta) + x_A
    y_mo[co_mul] = (radius + (randr2)) * np.sin(theta) + y_A
    xhe = np.hstack(x_mo[co_mul])
    yhe = np.hstack(y_mo[co_mul])
    i , j = 1 , 1
    x_cts[0] = xhe[0]
    y_cts[0] = yhe[0]
    #Aggregation coefficient for each step in heterogenous mode
    agg_co[0] = 1
    
    if agg > 0 :
        while light == 'red':
            prob_m = np.random.random()    
        
            if prob_m <= agg:
                pr_di = agg/i
                L[0] = agg_co[0]*pr_di #weight of every aggregation by all proteins
                if cts > 1:
                    for k  in range(cts - 1):
                        L[k + 1] = L[k] + agg_co[k + 1]*pr_di #weight of every aggregation by all proteins            
                
                    for co in range(cts):
                        if prob_m < L[co]:
                            ind = co
                            break
                direction = random.choice([-1, 1])        
                if direction == 1:
                    ran_p[ind] = ran_p[ind] + (d_pr / (2*np.pi*radius)) 
                    ropn = ran_p[ind]
                else:
                    ran_n[ind] = ran_n[ind] - (d_pr / (2*np.pi*radius))
                    ropn = ran_n[ind]
                xhe[i] = (radius + (sigma)) * np.cos(theta[ind] + ropn)  + x_A
                yhe[i] = (radius + (sigma)) * np.sin(theta[ind] + ropn)  + y_A
                agg_co[ind] = agg_co[ind] + 1
            else:
                x_cts[cts] = xhe[i]
                y_cts[cts] = yhe[i]
                agg_co[cts] = 1
                cts = cts + 1
        
            i = i + 1
        
            if i == Nr:
                light = 'green'  
    xhe = np.hstack(xhe)
    yhe = np.hstack(yhe)
    xt = []
    yt = []
    for i in range(int(Nr)):
        #Number of localizations
        Nr_fl = int(np.random.exponential(meann))
        if Nr_fl < 1:
            Nr_fl = 1
        x = np.array([None]*int(Nr_fl))
        y = np.array([None]*int(Nr_fl))
        mu_f, sigma_f = 0, 2*np.pi  
        

        randr1_f = np.array([None]*Nr_fl)
        randr2_f = np.array([None]*Nr_fl)


        randr1_f = np.random.normal(loc = 0, scale = sigma_fl, size = Nr_fl)
        randr2_f = np.random.normal(loc = 0, scale = sigma_fl, size = Nr_fl)
        x = (randr1_f) + xhe[i]
        y = (randr2_f) + yhe[i]    

        xt.append(x)
        yt.append(y)
    xf = np.hstack(xt) # final localisations
    yf = np.hstack(yt)
    x_wr[co_mul] = xf
    y_wr[co_mul] = yf 
    
    #Circle fit
    
    #Fitting to find the center in order to avoid cross section between droplets
    def calc_R(xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        aba =  np.sqrt((xhe-xc)**2 + (yhe-yc)**2)
        return aba

    def f_2(c):
        """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
        Ri = calc_R(*c)
        return Ri - Ri.mean()

    xm = np.mean(xhe, axis = 0)
    ym = np.mean(yhe, axis = 0)
    center_estimate = xm, ym
    center[co_mul], ier = optimize.leastsq(f_2, center_estimate)
    centerf[co_mul]= center[co_mul].reshape(1,2)
    xc_2[co_mul], yc_2[co_mul] = centerf[co_mul]
    Ri_fit       = calc_R(xc_2[co_mul],yc_2[co_mul])
    R_fit[co_mul]        = np.mean(Ri_fit)
    # Here we check if the droplets are interfere with each other or not
    if co_mul > 0:
        m = 1
        j = co_mul
        while m==1 :
            R_cr = np.sqrt((centerf[co_mul,0] - centerf[j - 1,0] )**2 + (centerf[co_mul,1] - centerf[j - 1,1])**2)
            if R_cr < R_fit[co_mul] + R_fit[j - 1] :
                co_mul = co_mul - 1
                m = 0;
            j = j - 1
            if j == 0:
              m = 0;  
    R_simu[co_mul] = radius
    co_mul = co_mul + 1
    
    if co_mul == NumDrop:
        flag = 'True'
#Data Extraction
# Ask for directory
root = tk.Tk()
root.withdraw()  # Hide the main window
infile = filedialog.askdirectory(title="Please choose your directory to save the simulated droplets data.")
if not infile:
    sys.exit("No directory selected or canceled by user.")

# Check if the directory exists
if not os.path.exists(infile):
    sys.exit("ERROR: The input directory does not exist.")
for n in range (NumDrop):
    col1 = "x [nm]"
    col2 = "y [nm]"
    data = pd.DataFrame({col1:x_wr[n],col2:y_wr[n]})
    nn = n + 1
    data.to_csv( infile + r'\Droplet_%d.csv'%nn,  index=False)

col1 = "R (nm)"
data = pd.DataFrame({col1:R_fit})
data.to_csv( infile + r'\R_fit.csv',  index=False)

# Plotting
plt.clf()
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
plt.savefig( infile + r'\FOV.png')
plt.show()


