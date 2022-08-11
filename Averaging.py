
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
#from sklearn import preprocessing

N = 2000

with open(r'C:\Users\JABER003\OneDrive - WageningenUR\Paper\Heterogeneity\Data\Simulation\Heterogeneous data\R_fit.csv') as file_name:
    R = np.loadtxt(file_name, delimiter=",", skiprows =1 )
NumDrop = len(R)
ch = np.array([None]*N)
ch[0] = 0;
Sum = 0
for i in range(N - 1):
    ch[i + 1] =  ch[i] + 0.0005
j = 0
for i in range(NumDrop):
    j = i + 1
    with open(r'C:\Users\JABER003\OneDrive - WageningenUR\Paper\Heterogeneity\Data\Simulation\Heterogeneous data\Relative_positions_%d.csv'%j) as file_name:
        dataf = np.loadtxt(file_name, delimiter=",", skiprows =1, usecols=2 )
    hh = pd.cut(dataf, bins = N)
    tt = hh.value_counts()
    v = tt.to_numpy()
    Sum = Sum + v
    Sum_max = np.amax(Sum)
    Sum2 = Sum/Sum_max
    print(j)
col1 = "ch"
col2 = "sum"
data = pd.DataFrame({col1:ch,col2:Sum2})
data.to_csv( r'C:\Users\JABER003\OneDrive - WageningenUR\Paper\Heterogeneity\Data\Simulation\Heterogeneous data\Avergae.csv',  index=False)
plt.plot(ch,Sum2)
plt.xlim([0,1])
plt.fill(ch, Sum2)
plt.savefig(r'C:\Users\JABER003\OneDrive - WageningenUR\Paper\Heterogeneity\Data\Simulation\Heterogeneous data\Average_RPD_normal.png')
plt.show()

