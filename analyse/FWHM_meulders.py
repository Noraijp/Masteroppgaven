import os
import numpy as np
import matplotlib.pyplot as plt

#def FWHM()
#path = '/Users/Nora/Documents/UiO/Masteroppgaven/'


file_33 = '/Users/Nora/Documents/UiO/Masteroppgaven/' + './Jon code/meulders_33MeV.csv'
file_16 = '/Users/Nora/Documents/UiO/Masteroppgaven/' + './Jon code/meulders_16MeV.csv'

with open(file_16) as f:
    begin = f.readlines()[1:]
    Energy_16MeV = []; N_counts_16MeV = []
    for line in begin:
        lines = line.split(',')
        Energy_16MeV.append(float(lines[0]))
        N_counts_16MeV.append(float(lines[1]))
    print('Energy_16MeV: ', Energy_16MeV)


with open(file_33) as p:
    begin2 = p.readlines()[1:]
    Energy_33MeV = []; N_counts_33MeV = []
    for line2 in begin2:
        lines2 = line2.split(',')
        Energy_33MeV.append(float(lines2[0]))
        N_counts_33MeV.append(float(lines2[1]))
