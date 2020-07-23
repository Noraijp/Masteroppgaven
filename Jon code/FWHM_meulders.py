import os
import numpy as np
import matplotlib.pyplot as plt

#def FWHM()
path_ = '/Users/Nora/Documents/UiO/Masteroppgaven/'

file_33 = os.getcwd() + '/meulders_33MeV.csv'
file_16 = os.getcwd() + '/meulders_16MeV.csv'

#file_33 = os.getcwd() + '/E33_flux_estimate.csv'
#file_16 = os.getcwd() + '/E16_Harrig_flux.csv'

with open(file_16) as f:
    begin = f.readlines()[1:]
    Energy_16MeV = []; N_counts_16MeV = []
    for line in begin:
        lines = line.split(',')
        Energy_16MeV.append(float(lines[0]))
        N_counts_16MeV.append(float(lines[1]))
    #print('Energy_16MeV: ', Energy_16MeV)
    #plt.plot(Energy_16MeV, N_counts_16MeV)
    #plt.show()


with open(file_33) as p:
    begin2 = p.readlines()[1:]
    Energy_33MeV = []; N_counts_33MeV = []
    for line2 in begin2:
        lines2 = line2.split(',')
        Energy_33MeV.append(float(lines2[0]))
        N_counts_33MeV.append(float(lines2[1]))
    max_y = max(N_counts_33MeV)
    #print('Maximum value of Y :', max_y)
    half_max_y = max_y / 2
    print('Halvpaten av max verdien for y: ', half_max_y)
    #max_x = Energy_33MeV[y_av.index(HM_y)]
    #plt.plot(Energy_33MeV, N_counts_33MeV)
    #x = np.interp(8550000000, Energy_33MeV, N_counts_33MeV)
    #print(x)
    #plt.plot(Energy_33MeV, N_counts_33MeV)
    #plt.show()


def lin_interp(x, y, i, half):
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

def half_max_F(x,y):
    half = max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return [lin_interp(x, y, zero_crossings_i[0], half),
                lin_interp(x, y, zero_crossings_i[1], half)]

def half_max_F_16(x,y):
    half = max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return [lin_interp(x, y, zero_crossings_i[0], half)]

E_33 = Energy_33MeV
F_33 = N_counts_33MeV

#print('!!___!!!!_!_!_!_!_!_!_!_')
#print(E_33)
#print(F_33)

E_16 = Energy_16MeV
F_16 = N_counts_16MeV

hmx_33 = half_max_F(E_33, F_33)
print('x1 and x2 for 33MeV: ', hmx_33)
#half_max.append(hmx)
#(mu,sigma) = norm.fit(E[i])
E_33 = np.array(E_33)
mu_33 = np.trapz(F_33*E_33, E_33)/np.trapz(F_33,E_33)
fwhm_33 = hmx_33[1]-hmx_33[0]
dEl_33 = mu_33-hmx_33[0]; dEr_33 = hmx_33[1]-mu_33   #left and right uncertainty in energy
#mu_array = mu
print(fwhm_33)
print('Right_33: ', dEr_33)
print('Left_33: ', dEl_33)
#print(type(hmx_33))


hmx_16 = half_max_F_16(E_16, F_16)
hmx_16.insert(0, 2.5311)
print('x1 and x2 for 16MeV: ', hmx_16)

print(type(hmx_16[0]))
print(type(hmx_16[1]))
#half_max.append(hmx)
#(mu,sigma) = norm.fit(E[i])
E_16 = np.array(E_16)
mu_16 = np.trapz(F_16*E_16, E_16)/np.trapz(F_16,E_16)
fwhm_16 = hmx_16[1]-hmx_16[0]
dEl_16 = mu_16-hmx_16[0]; dEr_16 = hmx_16[1]-mu_16   #left and right uncertainty in energy
#mu_array = mu
print(fwhm_16)
print('Right_16: ', dEr_16)
print('Left_16: ', dEl_16)
