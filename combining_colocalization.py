import os
import numpy as np
import pandas as pd


path = r'C:\Users\JABER003\OneDrive - WageningenUR\Paper\Heterogeneity\Data\EggYolk\Red_channel'
count = 0
object = os.scandir(path)
for n in object:
    count = count + 1
object.close()

for i in range(count):
    j = i + 1
    with open(r'C:\Users\JABER003\OneDrive - WageningenUR\Paper\Heterogeneity\Data\EggYolk\Red_channel\Droplet_%d.csv'%j) as file_name:
        datafR = np.loadtxt(file_name, delimiter=",", skiprows =1 )
    with open(r'C:\Users\JABER003\OneDrive - WageningenUR\Paper\Heterogeneity\Data\EggYolk\Green_channel\Droplet_%d.csv'%j) as file_name:
        datafG = np.loadtxt(file_name, delimiter=",", skiprows =1 )
    col1 = 'x_nm'
    col2 = 'y_nm'
    col3 = 'Error_nm'
    col4 = 'channel'     
    dataR = pd.DataFrame({col1:datafR[:,0],col2:datafR[:,1],col3:datafR[:,2],col4:1})
    dataG = pd.DataFrame({col1:datafG[:,0],col2:datafG[:,1],col3:datafG[:,2],col4:2})
    dataR.to_csv( r'C:\Users\JABER003\OneDrive - WageningenUR\Paper\Heterogeneity\Data\EggYolk\Combined_droplets\Droplet_%d.csv'%j,  index=False)
    dataG.to_csv( r'C:\Users\JABER003\OneDrive - WageningenUR\Paper\Heterogeneity\Data\EggYolk\Combined_droplets\Droplet_%d.csv'%j,  index=False, header=False, mode = "a")