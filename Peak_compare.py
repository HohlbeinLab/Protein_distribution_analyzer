import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

import pandas as pd
from tkinter import filedialog
import tkinter as tk


root = tk.Tk()
root.withdraw()  # Hide the main window
infile = filedialog.askdirectory(title="Please choose your directory to load the average histogram csv files.")
if not infile:
    sys.exit("No directory selected or canceled by user.")

if not os.path.exists(infile):
    sys.exit("ERROR; The input file does not exist.")

if infile[-4:] == '.csv':
    try:
        line = open(infile).readline()
        dataf = np.loadtxt(infile, delimiter=',',skiprows=1)
    except (EOFError, IOError, OSError) as exception:
        print("\n\nCould not open file: ", infile)
        print("\n\n", type(exception))
        sys.exit("Could not open the input file "+infile+".\n")
        
        
else:
    print('Sorry, wrong format!')
    sys.exit("The input file "+infile+" has the wrong format.")
x = dataf[:,1]
peaks,properties= find_peaks(x, prominence=(0, 10))
xp = [0,0]
nn = properties["prominences"]
mm = x[peaks[0]:peaks[-1]]
xp[0] = peaks[0]
xp[1] = peaks[-1]

#ymin= min(mm)
ymax = x[peaks]
plt.plot(x)
plt.plot(xp, x[xp], "x")
#plt.plot(peaks, x[peaks], "x")
path, in_file_no_path = os.path.split(infile)
plt.vlines(x=xp, ymin = 0,ymax = 1, color = "C1", linewidth=2.5, linestyle='dashed')
#plt.vlines(x=xp, ymin = ymin,ymax = x[xp], color = "C1")
#y_tot = ymax - ymin
#plt.hlines(y=ymin, xmin=peaks[0],xmax=peaks[-1], color = "C1")
plt.savefig(path + r"\Peaks_data.png")
plt.show()
col1 = "Distance"
col2 = "Peaks"
#data = pd.DataFrame({col1:xp,col2:x[xp] - ymin})
data = pd.DataFrame({col1:xp,col2:x[xp]})
data.to_csv( path + r"\Peaks_data2.csv",  index=False)
