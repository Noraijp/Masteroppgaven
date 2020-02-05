from math import *
sigma_A = 11.4/100
N_A = 6992
sigma_B = 1.1/100
N_B = 83329
sigma = sqrt(((sigma_A)*N_A)**2 + ((sigma_B)*N_B)**2)
N_tot = N_A + N_B
sigma_tot = sigma / N_tot
print(sigma)
print(N_tot)
print(sigma_tot)
