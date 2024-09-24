
import sys
import os
#Import PERPL to calculate the RPDs
sys.path.append(os.path.abspath(r'C:\Users\JABER003\OneDrive - WageningenUR\Programming\Heterogeneity analyzer\PERPL'))
from PERPL import relative_positions as rp
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from tkinter import filedialog
import tkinter as tk
#Bin size
N = 100

root = tk.Tk()
root.withdraw()  # Hide the main window
infile = filedialog.askdirectory(title="Please choose your directory to load the droplets' csv files.")
if not infile:
    sys.exit("No directory selected or canceled by user.")

if not os.path.exists(infile):
    sys.exit("ERROR; The input file does not exist.")


with open(infile + r'\R_fit.csv') as file_name:
        filterdist = np.loadtxt(file_name, delimiter=",", skiprows =1 )
        #Number of droplets
        NumDrop = len(filterdist)
        j = 0
        i = 0
#The diameter of each droplet        
Maxdist = [float(0) for x in range(NumDrop)]
  
for i in range(NumDrop):
    j = i + 1

    with open(infile + r'\Droplet_%d.csv'%j) as file_name:
        dataf = np.loadtxt(file_name, delimiter=",", skiprows =1 )
    dim = dataf.ndim
    #Get all distances between localizations
    nbr = rp.getdistances(dataf, filterdist[i]*3)
    rpd = rp.get_vectors(nbr, dim)
    # to be sure to collect all distances
    xx = np.where(rpd[:,3] > filterdist[i]*3)
    xx = np.hstack(xx)
    #rpd2 = rpd[0:,3]
    #Plotting the histogram to obtain the values for averagin
    if len(xx) != 0:
        ja  = plt.hist(rpd[0:xx[0],3] , bins = N)
    if len(xx) == 0:
        #xy_seperations that are plotted in histogram
        ja  = plt.hist(rpd[:,3] , bins = N)
    v = np.hstack(ja[0])
    Maxdist[i] = np.amax(rpd[:,3])
    #Bins channels
    ch = np.hstack(ja[1])
    ch.resize(N,)
    ch_max = ch[N - 1]
    ch = ch/ch_max *100
    #Normalized ja
    Sum_max = np.amax(v)
    Sum2 = v#/Sum_max
    #Exporting
    col1 = "ch"
    col2 = "sum"
    data = pd.DataFrame({col1:ch,col2:Sum2})
    data.to_csv(infile  + r'\Relative_positions_%d.csv'%j,  index=False)
    #Plotting
    plt.clf()
    plt.style.use('default')
    plt.plot(ch, Sum2)
    #plt.xlim([0,1])
    #plt.ylim([0,1])
    plt.savefig(infile + r'\Normalized_xy_%d.png'%j)
    plt.show()

# =============================================================================
# If you need all data frop RPD you can use these data
#         col1 = "xx_separation"
#         col2 = 'yy_separation'
#         col3 = 'xy_separation'
#         data = pd.DataFrame({col1:rpd[:,0],col2:rpd[:,1],col3:rpd[:,3]})
#         data.to_csv( r'C:\Users\JABER003\OneDrive - WageningenUR\Paper\Heterogeneity\Data\Simulation\Aggregation\%s\Relative_positions_%s.csv'% (ff[ab], j),  index=False)
# =============================================================================
    plt.clf()
    plt.style.use('default')
    plt.hist(rpd[:,3] , bins = N)
    plt.savefig(infile + r'\xy_separation_%d.png'%j)
    plt.clf()
    plt.style.use('default')
    plt.scatter(dataf[:,0], dataf[:,1], s = 5)
    plt.axis('equal')
    plt.savefig(infile +r'\Droplet_%d.png'%j)
    print(j)
data2 = pd.DataFrame({col1:Maxdist})
data2.to_csv(infile  + r'\Max_Distance.csv',  index=False)