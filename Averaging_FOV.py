
import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.special import erf
from tkinter.filedialog import askdirectory
from tkinter import filedialog
import tkinter as tk
#Number of Droplets
N = 100
root = tk.Tk()
root.withdraw()  # Hide the main window
infile = filedialog.askdirectory(title="Please choose your directory to load the relative position distributions of droplets (.csv)")
if not infile:
    sys.exit("No directory selected or canceled by user.")

if not os.path.exists(infile):
    sys.exit("ERROR; The input file does not exist.")
ch = [float(0) for x in range(N)]
ch = np.hstack(ch)
Sum = 0
j = 0
for i in range(N):
    j = i + 1
    with open(infile + r'\Relative_positions_%d.csv'%j) as file_name:
        dataf = np.loadtxt(file_name, delimiter=",", skiprows =1 )
    a = np.max(dataf[:,1])
    dataf[:,1] = dataf[:,1]/a
    Sum = Sum + dataf[:,1]
    print(j)
Sum_max = np.amax(Sum)
Sum2 = Sum/(Sum_max)


aaa = dataf[:,0]/100 
fig, ax = plt.subplots()
#Data extraction
col1 = "ch"
col2 = "sum"
data = pd.DataFrame({col1:dataf[:,0],col2:Sum2})
data.to_csv(infile + r'\Averge.csv',  index=False)
#Plotting
ax.plot(aaa,Sum2, fillstyle = 'full',color = 'cornflowerblue')
#ax.rcParams['font.size'] = 40
#ax.xlim([0,1])
#ax.ylim([0,1])
ax.set_xticks([0, 1])
ax.set_xticklabels(['0', '1'])
ax.set_yticks([0, 1])
ax.set_yticklabels(['0', '1'])
ax.fill_between(aaa, Sum2, color = 'cornflowerblue')
ax.tick_params(axis='both', which='major', labelsize=55)
#ax.tight_layout()
#plt.axis('off')
plt.gca().set_position([0, 0, 1, 1])

plt.savefig(infile + r'\Average_RPD_normal.png',bbox_inches='tight',transparent=True)
plt.show()

#FIT===============================================FIT===========================================================================FIT
