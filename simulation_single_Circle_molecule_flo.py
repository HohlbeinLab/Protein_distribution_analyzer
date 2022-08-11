# Python program to Plot Circle

# importing libraries
import numpy as np
from matplotlib import pyplot as plt
import scipy.stats as stats
from scipy import optimize
import random
# Creating equally spaced 100 data in range 0 to 2*pi


# Setting radius
m_density = 0.004 #1/nm
radius = 1000
sigma = radius/50
Nr = m_density*2*np.pi*radius
lower, upper = 0, 2*np.pi
mu_m, sigma_m = random.random()*2*np.pi, np.pi  
Mol = stats.truncnorm(
    (lower - mu_m) / sigma_m, (upper - mu_m) / sigma_m, loc=mu_m, scale=sigma_m)
theta = Mol.rvs(int(Nr))

# Generating x and y molecule position
x_mo = radius * np.cos(theta) + np.random.randn(int(Nr))*sigma
y_mo = radius * np.sin(theta) + np.random.randn(int(Nr))*sigma

F_density = 1/40
radius_m = 40
Nr_fl = F_density*2*np.pi*radius_m
mu_f, sigma_f = random.random()*2*np.pi, np.pi  
Fl = stats.truncnorm(
    (lower - mu_f) / sigma_f, (upper - mu_f) / sigma_f, loc=mu_f, scale=sigma_f)
sigma_fl = radius_m/5
theta_fl = Fl.rvs(int(Nr_fl))
x = np.array([None]*int(Nr))
y = np.array([None]*int(Nr))
for i in range(int(Nr)):
    x[i] = radius_m * np.cos(theta_fl) + np.random.randn(int(Nr_fl))*sigma_fl + x_mo[i]
    y[i] = radius_m * np.sin(theta_fl) + np.random.randn(int(Nr_fl))*sigma_fl + y_mo[i]
 

xf = np.hstack(x) # final localisations
yf = np.hstack(y)

#Circle fit

def calc_R(xc, yc):
     #calculate the distance of each 2D points from the center (xc, yc) 
    return np.sqrt((xf-xc)**2 + (yf-yc)**2)

def f_2(c):
     #calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) 
    Ri = calc_R(*c)
    return Ri - Ri.mean()

xm = np.mean(xf, axis = 0)
ym = np.mean(yf, axis = 0)
center_estimate = xm, ym
center_2, ier = optimize.leastsq(f_2, center_estimate)
xc_2, yc_2 = center_2
Ri_2       = calc_R(*center_2)
R_2        = np.mean(Ri_2)

print(center_2)

# Plotting
figure, axs = plt.subplots() 
#plt.scatter(x_mo, y_mo , s = 10, alpha = 0.5, zorder = 1) # plot molecules position
plt.scatter(xf, yf , s = 10, alpha = 0.5, zorder = 1) # plot fluorophore
''' # color coded version
for i in range(int(Nr)):
    if i%2 == 0:
        plt.scatter(x[i], y[i] , s = 10, alpha = 0.5, zorder = 1, color = 'b')
    else:
        plt.scatter(x[i], y[i] , s = 10, alpha = 0.5, zorder = 1, color = 'r')
'''
cc = plt.Circle(center_2, R_2, fill=False, linestyle='--', linewidth = 2, color = 'tab:orange', zorder = 2) 


axs.add_patch( cc ) 
axs.text(0.98, 0.005, f"R = %d nm"%R_2 ,
        verticalalignment='bottom', horizontalalignment='right',
        transform=axs.transAxes,
        color='tab:orange', fontsize=10)



plt.axis('equal')
plt.title('Circle')

plt.show()