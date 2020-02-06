import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#filename = "64Cu_230.dat";lambda_ = np.log(2)/(12.701)#*60)
filename = "24Na_Al16MeV_113.dat";lambda_ = np.log(2)/(14.997)#*60)

time = np.genfromtxt(filename, delimiter = ",", usecols = [0])
A = np.genfromtxt(filename, delimiter = ",", usecols = [1])
sigma_A = np.genfromtxt(filename, delimiter = ",", usecols = [2])

index = ~(np.isnan(A) | np.isnan(sigma_A))


def direct_decay(time, A0_guess):
    A_est = A0_guess*np.exp(-lambda_*time)
    return A_est

popt, pcov = curve_fit(direct_decay, time[index], A[index], p0=600, sigma=sigma_A[index], absolute_sigma=True)
print(popt)
xplot = np.linspace(0,np.max(time), 1000)

plt.plot(time[index], A[index],".")
plt.plot(xplot, direct_decay(xplot, *popt))
plt.show()
