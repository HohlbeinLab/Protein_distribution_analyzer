
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from sklearn import preprocessing


with open(r'C:\Users\JABER003\OneDrive - WageningenUR\Programming\Heterogeneity analyzer\Data_simulation\Heterogeneous data\Relative_positions_75.csv') as file_name:
    dataf = np.loadtxt(file_name, delimiter=",", skiprows =1, usecols=2 )
