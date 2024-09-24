import scipy.special as ss
import numpy as np
from matplotlib import pyplot as plt
import sys
import os
sys.path.append(os.path.abspath(r'C:\Users\JABER003\OneDrive - WageningenUR\Programming\Heterogeneity analyzer\PERPL'))
from PERPL import relative_positions as rp
from relative_positions import getdistances
from relative_positions import get_vectors
from scipy import i0
#Diameter of a droplet
D = 1000
#Number of proteins
Nr = 200
##Relative broadening
sigma =60

v = np.zeros((Nr, 2))

su = [float(0) for i in range(Nr)]
su = np.hstack(su)
sumn = np.hstack(su)
ch = [float(0) for i in range(Nr)]
ch = np.hstack(ch)
vf = np.hstack(ch)
j = 0
flag = 'b'
for i in range(Nr):
    #define a structure of the droplet (ring shape)
    v[i,0] = (D/2) * np.cos(2*np.pi*i/Nr)
    v[i,1] = (D/2) * np.sin(2*np.pi*i/Nr)
#calculating the relative positions
relative_positions = getdistances(v, D, sort_and_halve=False)
xy_separations = np.sqrt(relative_positions[:, 0] ** 2
                             + relative_positions[:, 1] ** 2)
#The defined function for non gausian distribution

def func2D(r, rmean, sigma):
    
    if np.isscalar(r):
        p = 0 # Might help to diagnose if there are problems
        if (rmean * r / sigma ** 2) < 700:   
            p = (r / sigma ** 2) * (np.exp(-(rmean ** 2 + r ** 2) \
                / (2 * sigma ** 2)) * i0(r * rmean / sigma ** 2))
        else:  # Approximate overly large i0()
            p = 1 / (np.sqrt(2 * np.pi) * sigma) * np.sqrt(r / rmean) \
                * np.exp(-((r - rmean) ** 2) / (2 * sigma ** 2))

    # If z is an array of distances at which to evaluate the function:    
    else:
        if (np.max(r) * rmean / sigma ** 2) < 700.:
            p = (r / sigma ** 2
                    * (np.exp(-(rmean ** 2 + r ** 2)
                                / (2 * sigma ** 2)
                                )
                        * i0(r * rmean / sigma ** 2)
                        )
                    )
        else:  # Approximate for overly large i0()
            p = (1 / (np.sqrt(2 * np.pi) * sigma)
                    * np.sqrt(r / rmean)
                    * np.exp(-((r - rmean) ** 2) / (2 * sigma ** 2))
                    )
    return p
j = 0
while flag == 'b':
    su = [float(0) for i in range(Nr)]
    su = np.hstack(su)    
    
                       

    for i in range(Nr):
        #for k in range (Nr):
            su[i] =su[i] + func2D(xy_separations[i], xy_separations[j] , sigma) 
    
    sumn = sumn + su
    j = j + 1
    if j == Nr:
        flag = 'a'
#plt.hist(xy_separations)        
for i in range(Nr):
    ch[i] = xy_separations[i]
#plt.hist(su, bins = 100)
plt.plot(ch, sumn)
plt.show()



















