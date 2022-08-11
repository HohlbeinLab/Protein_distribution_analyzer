# -*- coding: utf-8 -*-
"""
Created on Wed May  4 19:05:06 2022

@author: jaber003
"""
import sys
import os
sys.path.append(os.path.abspath(r'C:\Users\JABER003\OneDrive - WageningenUR\Programming\Heterogeneity analyzer\PERPL'))
from PERPL import relative_positions as rp
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd



with open(r'C:\Users\JABER003\OneDrive - WageningenUR\Programming\Heterogeneity analyzer\Data_simulation\R_fit.csv') as file_name:
    filterdist = np.loadtxt(file_name, delimiter=",", skiprows =1 )
NumDrop = len(filterdist)
j = 0
i = 0

for i in range(NumDrop):
    j = i + 1
    with open(r'C:\Users\JABER003\OneDrive - WageningenUR\Programming\Heterogeneity analyzer\Data_simulation\Homogeneous data (low unc)\Droplet_%d.csv'%j) as file_name:
        dataf = np.loadtxt(file_name, delimiter=",", skiprows =1 )
    dim = dataf.ndim
    nbr = rp.getdistances(dataf, filterdist[i]*2.1)
    rpd = rp.get_vectors(nbr, dim)
    col1 = "xx_separation"
    col2 = 'yy_separation'
    col3 = 'xy_separation'
    data = pd.DataFrame({col1:rpd[:,0],col2:rpd[:,1],col3:rpd[:,3]})
    data.to_csv( r'C:\Users\JABER003\OneDrive - WageningenUR\Programming\Heterogeneity analyzer\Data_simulation\Homogeneous data (low unc)\Relative_positions_%d.csv'%j,  index=False)
    plt.clf()
    plt.style.use('ggplot')
    plt.hist(rpd[:,3] , bins = 2000)
    plt.savefig(r'C:\Users\JABER003\OneDrive - WageningenUR\Programming\Heterogeneity analyzer\Data_simulation\Homogeneous data (low unc)\xy_separation_%d.png'%j)
    plt.clf()
    plt.style.use('default')
    plt.scatter(dataf[:,0], dataf[:,1], s = 1)
    plt.axis('equal')
    plt.savefig(r'C:\Users\JABER003\OneDrive - WageningenUR\Programming\Heterogeneity analyzer\Data_simulation\Homogeneous data (low unc)\Droplet_%d.png'%j)
    print(j)