import numpy as np
import matplotlib.pyplot as plt
import npat
import csv

from npat import *


with open('meulders.csv') as f:
	meulders_33MeV = np.array([i.split(',') for i in f.read().split('\n')[:-1]], dtype=np.float64)

# plt.plot(meulders_33MeV[:,0],meulders_33MeV[:,1])
# plt.show()


### Read in numbers of decays from csv file
def decomment(csvfile):
    for row in csvfile:
        raw = row.split('#')[0].strip()
        if raw: yield raw



def read_csv(name_of_csv_file):
	results = []
	with open('../analyse/find activity/activity_csv/'+name_of_csv_file) as csvfile:
	    reader = csv.reader(decomment(csvfile))
	    for row in reader:
	        results.append(row)

	return np.asarray(results, dtype=float)


def calculate_flux(csv_list, reaction_list, target_list, product_list, mass_33MeV, unc_mass_33MeV):
	flux_array  =  []
	unc_flux_array  =  []

	for i in range(len(reaction_list)):
		activity_data = read_csv(csv_list[i])
		#print('CSV data: \n',activity_data)

		activity_16MeV = activity_data[:,0]
		activity_33MeV = activity_data[:,1]
		#print('16 MeV data: \n',activity_16MeV)
		#print('33 MeV data: \n',activity_33MeV)


		### starting with yttrium monitor
		## Only for monitor reactions

		rx_33MeV = npat.Reaction(reaction_list[i],library='IRDFF')
		# rx_33MeV.plot(scale='loglog')

		xs_avg_33MeV, unc_xs_avg_33MeV = rx_33MeV.average(meulders_33MeV[:,0], meulders_33MeV[:,1], unc=True)
		#xs_avg_33MeV *= 1E-3  ## convert mb to b
		print('Meulders 33 MeV averaged XS (mb):', xs_avg_33MeV)



		dc_33MeV =      DecayChain(product_list[i], R=1, time=530.0)
		dc_33MeV.append(DecayChain(product_list[i],      time=((5.0*60.0)+40.0)))
		dc_33MeV.append(DecayChain(product_list[i], R=1, time=6750.0))

#		dc_33MeV = npat.DecayChain(product_list[i], R=1, time=irr_time_33MeV)
		dc_33MeV.append(npat.DecayChain(product_list[i], time=1E6))

		count_start_time = 0.0 # h
		measured_monitor_activity_33MeV = activity_33MeV[0]

		R_33MeV = measured_monitor_activity_33MeV/dc_33MeV.activity(product_list[i], t=count_start_time, units='h')
		unc_R_33MeV = activity_33MeV[1]/dc_33MeV.activity(product_list[i], t=count_start_time, units='h')



		target_isotope = npat.Isotope(target_list[i])
		ab = 1E-2*target_isotope.abundance()
		n_atoms_33MeV = (ab*mass_33MeV*6.022E-1)/target_isotope.mass

		avg_flux_33MeV = R_33MeV/(n_atoms_33MeV*xs_avg_33MeV)
		#unc_avg_flux_33MeV = np.sqrt(  (unc_R_33MeV/R_33MeV)**2  + (unc_xs_avg_33MeV/xs_avg_33MeV)**2 + (unc_mass_33MeV/mass_33MeV)**2   )
		unc_avg_flux_33MeV = avg_flux_33MeV*np.sqrt(  (activity_33MeV[1]/activity_33MeV[0])**2  + (unc_xs_avg_33MeV/xs_avg_33MeV)**2 + (unc_mass_33MeV/mass_33MeV)**2   )
		# print('% uncertainty in R: ', 100*activity_33MeV[1]/activity_33MeV[0])
		# print('% uncertainty in N_atoms: ', 100*unc_mass_33MeV/mass_33MeV)
		# print('% uncertainty in xs_avg: ', 100*unc_xs_avg_33MeV/xs_avg_33MeV)

		print('\n************************************************\n')
		print('33MeV average flux for monitor ',  reaction_list[i], ':  ', avg_flux_33MeV, ' +/- ', unc_avg_flux_33MeV, ' (', 100*unc_avg_flux_33MeV/avg_flux_33MeV, ' %)')
		flux_array.append(avg_flux_33MeV)
		unc_flux_array.append(unc_avg_flux_33MeV)

	print('\n************************************************\n')
	return flux_array, unc_flux_array






irr_time_33MeV = 1.5*3600.0 #sec


### for yttrium

csv_list      = ['Y_88Y.csv']
reaction_list = ['89Y(n,2n)88Y']
target_list   = ['89Y']
product_list  = ['88Y']
mass_33MeV = 0.4657 #g
unc_mass_33MeV = 0.0006  #g

y_avg_flux, y_unc_avg_flux = calculate_flux(csv_list, reaction_list, target_list, product_list, mass_33MeV, unc_mass_33MeV)
print('Y flux            : ', y_avg_flux)
print('Y flux uncertinty : ', y_unc_avg_flux)



### for indium

csv_list      = ['In_114mIn.csv', 'In_113mIn.csv',  'In_115mIn.csv',  'In_116mIn.csv']
reaction_list = ['113IN(n,g)114INm', '113IN(n,inl)113INm', '115IN(n,inl)115INm', '115IN(n,g)116INm']
target_list   = ['113IN', '113IN', '115IN', '115IN']
product_list  = ['114INm', '113INm', '115INm', '116INm']
mass_33MeV = 0.5443 #g
unc_mass_33MeV = 0.0013  #g

### Search the IRDFF neutron library for reactions
#lb = Library('IRDFF')
# print(lb.search(target='113IN'))
#print(lb.search(product='116INm'))
# print(lb.search(target='113IN',product='114INm1'))

in_avg_flux, in_unc_avg_flux = calculate_flux(csv_list, reaction_list, target_list, product_list, mass_33MeV, unc_mass_33MeV)
print('In flux            : ', in_avg_flux)
print('In flux uncertinty : ', in_unc_avg_flux)




### for Aluminum

csv_list      = ['Al_24Na.csv']
reaction_list = ['27AL(n,a)24NA']
target_list   = ['27AL']
product_list  = ['24NA']
mass_33MeV = 0.2563 #g
unc_mass_33MeV = 0.0015  #g

### Search the IRDFF neutron library for reactions
lb = Library('IRDFF')
# print(lb.search(target='113IN'))
print(lb.search(product='24Na'))
# print(lb.search(target='113IN',product='114INm1'))

al_avg_flux, al_unc_avg_flux = calculate_flux(csv_list, reaction_list, target_list, product_list, mass_33MeV, unc_mass_33MeV)
print('Al flux            : ', al_avg_flux)
print('Al flux uncertinty : ', al_unc_avg_flux)




### for Zirconium

csv_list      = ['Zr_89Zr.csv']
reaction_list = ['90ZR(n,2n)89ZR']
target_list   = ['90ZR']
product_list  = ['89ZR']
mass_33MeV = 0.7557 #g
unc_mass_33MeV = 0.0012  #g

### Search the IRDFF neutron library for reactions
#lb = Library('IRDFF')
# print(lb.search(target='113IN'))
#print(lb.search(product='89Zr'))
# print(lb.search(target='113IN',product='114INm1'))

zr_avg_flux, zr_unc_avg_flux = calculate_flux(csv_list, reaction_list, target_list, product_list, mass_33MeV, unc_mass_33MeV)
print('Zr flux            : ', zr_avg_flux)
print('Zr flux uncertinty : ', zr_unc_avg_flux)


all_fluxes = []
all_unc_fluxes = []

for i in (y_avg_flux, in_avg_flux, al_avg_flux, zr_avg_flux):
	all_fluxes.extend(i[:])
for i in (y_unc_avg_flux, in_unc_avg_flux, al_unc_avg_flux, zr_unc_avg_flux):
	all_unc_fluxes.extend(i[:])
print(all_fluxes)
# print(all_fluxes[3])
x_index = range(len(all_fluxes))


approximate_average_flux = np.average(all_fluxes, weights=1/np.square(all_unc_fluxes))
unc_approximate_average_flux = np.average(all_unc_fluxes)
print('approximate_average_flux: ',approximate_average_flux)

# plt.errorbar(x_index, all_fluxes, yerr=all_unc_fluxes, marker='.', linestyle='')
# plt.plot(x_index, approximate_average_flux*np.ones(len(x_index)), color='red')
#plt.show()





### for zinc
product_name='67CU'

product_activity_data = read_csv('Zn_67Cu.csv')
#print('CSV data: \n',activity_data)

product_activity_16MeV = product_activity_data[:,0]
product_activity_33MeV = product_activity_data[:,1]


product_dc_33MeV =      DecayChain(product_name, R=1, time=530.0)
product_dc_33MeV.append(DecayChain(product_name,      time=((5.0*60.0)+40.0)))
product_dc_33MeV.append(DecayChain(product_name, R=1, time=6750.0))

#dc = npat.DecayChain(product_name, R=1, time=irr_time)
product_dc_33MeV.append(npat.DecayChain(product_name, time=1E6))

count_start_time = 0.0 # h
#measured_A0_67ZN_activity = 3.7E3 # bq

#R = product_activity_33MeV[0]/product_dc_33MeV.activity(product_name, t=count_start_time, units='h')
R = product_activity_33MeV[0]

mass = 0.8427 #g
itp = npat.Isotope('67ZN')
ab = 1E-2*itp.abundance()

n_atoms = (ab*mass*6.022E-1)/itp.mass


xs = R/(n_atoms*approximate_average_flux)
print('Our measured cross section (mb); ', xs)

lb = Library('IRDFF')
#print(lb.search(target='113IN'))
print(lb.search(product='67CU'))
rx = npat.Reaction('67ZN(n,p)67CU', library='IRDFF')
meulders_33MeV_xs_avg = rx.average(meulders_33MeV[:,0], meulders_33MeV[:,1])
print('meulders_33MeV_xs_avg: ',meulders_33MeV_xs_avg)
