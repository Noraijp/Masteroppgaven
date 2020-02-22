import numpy as np
import matplotlib.pyplot as plt
import npat
import csv

from npat import *

from weighted_average import *

with open('meulders_33MeV.csv') as f:
	meulders_33MeV = np.array([i.split(',') for i in f.read().split('\n')[:-1]], dtype=np.float64)

with open('meulders_16MeV.csv') as f:
	meulders_16MeV = np.array([i.split(',') for i in f.read().split('\n')[:-1]], dtype=np.float64)

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


def calculate_flux(csv_list, reaction_list, target_list, product_list, mass_16MeV, unc_mass_16MeV, mass_33MeV, unc_mass_33MeV):
	flux_array_33MeV  =  []
	unc_flux_array_33MeV  =  []
	flux_array_16MeV  =  []
	unc_flux_array_16MeV  =  []
	production_rate =  np.zeros((2,len(reaction_list)))
	number_of_atoms =  np.zeros((2,len(reaction_list)))
	flux_avg_cross_section =  np.zeros((2,len(reaction_list)))
	unc_production_rate =  np.zeros((2,len(reaction_list)))
	unc_number_of_atoms =  np.zeros((2,len(reaction_list)))
	unc_flux_avg_cross_section =  np.zeros((2,len(reaction_list)))


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
		rx_16MeV = npat.Reaction(reaction_list[i],library='IRDFF')

		# rx_33MeV.plot(scale='loglog')

		xs_avg_33MeV, unc_xs_avg_33MeV = rx_33MeV.average(meulders_33MeV[:,0], meulders_33MeV[:,1], unc=True)
		#xs_avg_33MeV *= 1E-3  ## convert mb to b
		# print('Meulders 33 MeV averaged XS (mb):', xs_avg_33MeV)

		xs_avg_16MeV, unc_xs_avg_16MeV = rx_16MeV.average(meulders_16MeV[:,0], meulders_16MeV[:,1], unc=True)
		#xs_avg_33MeV *= 1E-3  ## convert mb to b
		# print('Meulders 16 MeV averaged XS (mb):', xs_avg_16MeV)


		dc_33MeV =      DecayChain(product_list[i], R=1, time=530.0)
		dc_33MeV.append(DecayChain(product_list[i],      time=((5.0*60.0)+40.0)))
		dc_33MeV.append(DecayChain(product_list[i], R=1, time=6750.0))


		dc_16MeV =      DecayChain(product_list[i], R=1, time= (7*3600) + (57*60))

#		dc_33MeV = npat.DecayChain(product_list[i], R=1, time=irr_time_33MeV)
		dc_33MeV.append(npat.DecayChain(product_list[i], time=1E6))
		dc_16MeV.append(npat.DecayChain(product_list[i], time=1E6))


		count_start_time = 0.0 # h
		measured_monitor_activity_33MeV = activity_33MeV[0]
		measured_monitor_activity_16MeV = activity_16MeV[0]

		R_33MeV = measured_monitor_activity_33MeV/dc_33MeV.activity(product_list[i], t=count_start_time, units='h')
		unc_R_33MeV = activity_33MeV[1]/dc_33MeV.activity(product_list[i], t=count_start_time, units='h')

		R_16MeV = measured_monitor_activity_16MeV/dc_16MeV.activity(product_list[i], t=count_start_time, units='h')
		unc_R_16MeV = activity_16MeV[1]/dc_16MeV.activity(product_list[i], t=count_start_time, units='h')



		target_isotope = npat.Isotope(target_list[i])
		#print('target_isotope: ', str(target_isotope)[0:3])
		if str(target_isotope)[0:3] == 'nat':
			ab = 1.0
		else:
			ab = 1E-2*target_isotope.abundance()
		#print('ab: ', ab)
		#print(target_isotope)
		if str(target_isotope) == 'natIN':
			isotope_mass = 114.818
		elif str(target_isotope) == 'natY':
			isotope_mass = 88.90584
		elif str(target_isotope) == 'natAL':
			isotope_mass = 26.9815384
		elif str(target_isotope) == 'natZR':
			isotope_mass = 91.224
		elif str(target_isotope) == 'natZN':
			isotope_mass = 65.38
		else:
			isotope_mass = target_isotope.mass
		n_atoms_33MeV = (ab*mass_33MeV*6.022E-1)/isotope_mass
		#print('n_atoms: ', n_atoms_33MeV)
		#print('target_mass: ', isotope_mass)
		n_atoms_16MeV = (ab*mass_16MeV*6.022E-1)/isotope_mass

		avg_flux_33MeV = R_33MeV/(n_atoms_33MeV*xs_avg_33MeV)
		avg_flux_16MeV = R_16MeV/(n_atoms_16MeV*xs_avg_16MeV)
		#unc_avg_flux_33MeV = np.sqrt(  (unc_R_33MeV/R_33MeV)**2  + (unc_xs_avg_33MeV/xs_avg_33MeV)**2 + (unc_mass_33MeV/mass_33MeV)**2   )
		unc_avg_flux_33MeV = avg_flux_33MeV*np.sqrt(  (activity_33MeV[1]/activity_33MeV[0])**2  + (unc_xs_avg_33MeV/xs_avg_33MeV)**2 + (unc_mass_33MeV/mass_33MeV)**2   )
		unc_avg_flux_16MeV = avg_flux_16MeV*np.sqrt(  (activity_16MeV[1]/activity_16MeV[0])**2  + (unc_xs_avg_16MeV/xs_avg_16MeV)**2 + (unc_mass_16MeV/mass_16MeV)**2   )

		# print('% uncertainty in R: ', 100*activity_33MeV[1]/activity_33MeV[0])
		# print('% uncertainty in N_atoms: ', 100*unc_mass_33MeV/mass_33MeV)
		# print('% uncertainty in xs_avg: ', 100*unc_xs_avg_33MeV/xs_avg_33MeV)


		# print('\n************************************************\n')
		# print('33MeV average flux for monitor ',  reaction_list[i], ':  ', avg_flux_33MeV, ' +/- ', unc_avg_flux_33MeV, ' (', 100*unc_avg_flux_33MeV/avg_flux_33MeV, ' %)')
		flux_array_33MeV.append(avg_flux_33MeV)
		unc_flux_array_33MeV.append(unc_avg_flux_33MeV)

		# print('16MeV average flux for monitor ',  reaction_list[i], ':  ', avg_flux_16MeV, ' +/- ', unc_avg_flux_16MeV, ' (', 100*unc_avg_flux_16MeV/avg_flux_16MeV, ' %)')
		flux_array_16MeV.append(avg_flux_16MeV)
		unc_flux_array_16MeV.append(unc_avg_flux_16MeV)

		production_rate[:,i] = np.array((R_16MeV, R_33MeV))
		number_of_atoms[:,i] = np.array((n_atoms_16MeV, n_atoms_33MeV))
		flux_avg_cross_section[:,i] = np.array((xs_avg_16MeV, xs_avg_33MeV))
		unc_production_rate[:,i] = np.array((unc_R_16MeV, unc_R_33MeV))
		unc_number_of_atoms[:,i] = np.array((n_atoms_16MeV*(unc_mass_16MeV/mass_16MeV), n_atoms_33MeV*(unc_mass_33MeV/mass_33MeV)))
		unc_flux_avg_cross_section[:,i] = np.array((unc_xs_avg_16MeV, unc_xs_avg_33MeV))

	# print('\n************************************************\n')
	return flux_array_33MeV, unc_flux_array_33MeV, flux_array_16MeV, unc_flux_array_16MeV, production_rate, number_of_atoms, flux_avg_cross_section, unc_production_rate, unc_number_of_atoms, unc_flux_avg_cross_section


# def Average_BeamCurrent(production_rate, number_of_atoms, flux_avg_cross_section, unc_production_rate, unc_number_of_atoms, unc_flux_avg_cross_section, csv_filename='averaged_currents.csv'):






## for yttrium

csv_list      = ['Y_88Y.csv']
reaction_list = ['89Y(n,2n)88Y']
target_list   = ['89Y']
product_list  = ['88Y']
mass_33MeV = 0.4657 #g
mass_16MeV = 0.5053
unc_mass_33MeV = 0.0006  #g
unc_mass_16MeV = 0.0021

y_avg_flux_33MeV, y_unc_avg_flux_33MeV, y_avg_flux_16MeV, y_unc_avg_flux_16MeV, y_production_rate, y_number_of_atoms, y_flux_avg_cross_section, y_unc_production_rate, y_unc_number_of_atoms, y_unc_flux_avg_cross_section = calculate_flux(csv_list, reaction_list, target_list, product_list, mass_16MeV, unc_mass_16MeV, mass_33MeV, unc_mass_33MeV)
# print('Y flux at 33MeV           : ', y_avg_flux_33MeV)
# print('Y flux uncertinty_33MeV : ', y_unc_avg_flux_33MeV)
# print('Y flux at 16MeV           : ', y_avg_flux_16MeV)
# print('Y flux uncertinty_16MeV : ', y_unc_avg_flux_16MeV)



### for indium

csv_list      = ['In_114mIn.csv','In_114mIn.csv', 'In_113mIn.csv',  'In_115mIn.csv']
#reaction_list = ['113IN(n,g)114INm', '113IN(n,inl)113INm', '115IN(n,inl)115INm', '115IN(n,g)116INm']
reaction_list = ['natIN(n,x)114INm',  '115IN(n,2n)114INm', '113IN(n,inl)113INm', '115IN(n,inl)115INm']
target_list   = ['natIN',  '115IN', '113IN', '115IN']
product_list  = ['114INm',  '114INm','113INm', '115INm']
# csv_list      = ['In_114mIn.csv','In_114mIn.csv' ,'In_114mIn.csv', 'In_113mIn.csv',  'In_115mIn.csv',  'In_116mIn.csv']
# #reaction_list = ['113IN(n,g)114INm', '113IN(n,inl)113INm', '115IN(n,inl)115INm', '115IN(n,g)116INm']
# reaction_list = ['natIN(n,x)114INm', '113IN(n,g)114INm', '115IN(n,2n)114INm', '113IN(n,inl)113INm', '115IN(n,inl)115INm', '115IN(n,g)116INm']
# target_list   = ['natIN', '113IN', '115IN', '113IN', '115IN', '115IN']
# product_list  = ['114INm', '114INm', '114INm','113INm', '115INm', '116INm']
mass_33MeV = 0.5443 #g
unc_mass_33MeV = 0.0013  #g
mass_16MeV = 0.5530 #g
unc_mass_16MeV = 0.0035  #g

# ## Search the IRDFF neutron library for reactions
# lb = Library('IRDFF')
# print(lb.search(target='113IN'))
# print(lb.search(product='114INm'))
# #print(lb.search(target='113IN',product='114INm1'))

in_avg_flux_33MeV, in_unc_avg_flux_33MeV, in_avg_flux_16MeV, in_unc_avg_flux_16MeV, in_production_rate, in_number_of_atoms, in_flux_avg_cross_section, in_unc_production_rate, in_unc_number_of_atoms, in_unc_flux_avg_cross_section = calculate_flux(csv_list, reaction_list, target_list, product_list, mass_16MeV, unc_mass_16MeV, mass_33MeV, unc_mass_33MeV)
# print('In flux at 33MeV           : ', in_avg_flux_33MeV)
# print('In flux uncertinty_33MeV : ', in_unc_avg_flux_33MeV)
# print('In flux at 16MeV           : ', in_avg_flux_16MeV)
# print('In flux uncertinty_16MeV : ', in_unc_avg_flux_16MeV)




### for Aluminum

csv_list      = ['Al_24Na.csv', 'Al_24Na.csv']
reaction_list = ['27AL(n,x)24NA', '27AL(n,a)24NA']
target_list   = ['27AL', '27AL']
product_list  = ['24NA', '24NA']
mass_33MeV = 0.2563 #g
unc_mass_33MeV = 0.0015  #g
mass_16MeV = 0.2573 #g
unc_mass_16MeV = 0.0006 #g

# ### Search the IRDFF neutron library for reactions
# lb = Library('IRDFF')
# # print(lb.search(target='113IN'))
# print(lb.search(product='24Na'))
# # print(lb.search(target='113IN',product='114INm1'))

al_avg_flux_33MeV, al_unc_avg_flux_33MeV, al_avg_flux_16MeV, al_unc_avg_flux_16MeV, al_production_rate, al_number_of_atoms, al_flux_avg_cross_section, al_unc_production_rate, al_unc_number_of_atoms, al_unc_flux_avg_cross_section = calculate_flux(csv_list, reaction_list, target_list, product_list, mass_16MeV, unc_mass_16MeV, mass_33MeV, unc_mass_33MeV)
# print('Al flux at 33MeV           : ', al_avg_flux_33MeV)
# print('Al flux uncertinty_33MeV : ', al_unc_avg_flux_33MeV)
# print('Al flux at 16MeV           : ', al_avg_flux_16MeV)
# print('Al flux uncertinty_16MeV : ', al_unc_avg_flux_16MeV)




### for Zirconium

csv_list      = ['Zr_89Zr.csv', 'Zr_89Zr.csv']
reaction_list = ['natZR(n,x)89ZR', '90ZR(n,2n)89ZR']
target_list   = ['natZR','90ZR']
product_list  = ['89ZR', '89ZR']
mass_33MeV = 0.7557 #g
unc_mass_33MeV = 0.0012  #g
mass_16MeV = 0.7560 #g
unc_mass_16MeV = 0.0010 #g

# ## Search the IRDFF neutron library for reactions
# lb = Library('IRDFF')
# # print(lb.search(target='113IN'))
# print(lb.search(product='89Zr'))
# # print(lb.search(target='113IN',product='114INm1'))

zr_avg_flux_33MeV, zr_unc_avg_flux_33MeV, zr_avg_flux_16MeV, zr_unc_avg_flux_16MeV, zr_production_rate, zr_number_of_atoms, zr_flux_avg_cross_section, zr_unc_production_rate, zr_unc_number_of_atoms, zr_unc_flux_avg_cross_section = calculate_flux(csv_list, reaction_list, target_list, product_list, mass_16MeV, unc_mass_16MeV, mass_33MeV, unc_mass_33MeV)
# print('Zr flux at 33MeV           : ', zr_avg_flux_33MeV)
# print('Zr flux uncertinty_33MeV : ', zr_unc_avg_flux_33MeV)
# print('Zr flux at 16MeV           : ', zr_avg_flux_16MeV)
# print('Zr flux uncertinty_16MeV : ', zr_unc_avg_flux_16MeV)
#

all_fluxes_33MeV = []
all_unc_fluxes_33MeV = []

for i in (y_avg_flux_33MeV, in_avg_flux_33MeV, al_avg_flux_33MeV, zr_avg_flux_33MeV):
	all_fluxes_33MeV.extend(i[:])
for i in (y_unc_avg_flux_33MeV, in_unc_avg_flux_33MeV, al_unc_avg_flux_33MeV, zr_unc_avg_flux_33MeV):
	all_unc_fluxes_33MeV.extend(i[:])
# print(all_fluxes_33MeV)
# print(all_fluxes[3])
x_index_33MeV = range(len(all_fluxes_33MeV))

approximate_average_flux_33MeV = np.average(all_fluxes_33MeV, weights=1/np.square(all_unc_fluxes_33MeV))
unc_approximate_average_flux_33MeV = np.average(all_unc_fluxes_33MeV)
# print('approximate_average_flux_33MeV: ',approximate_average_flux_33MeV)


all_fluxes_16MeV = []
all_unc_fluxes_16MeV = []

for i in (y_avg_flux_16MeV, in_avg_flux_16MeV, al_avg_flux_16MeV, zr_avg_flux_16MeV):
	all_fluxes_16MeV.extend(i[:])
for i in (y_unc_avg_flux_16MeV, in_unc_avg_flux_16MeV, al_unc_avg_flux_16MeV, zr_unc_avg_flux_16MeV):
	all_unc_fluxes_16MeV.extend(i[:])
# print(all_fluxes_16MeV)
# print(all_fluxes[3])
x_index_16MeV = range(len(all_fluxes_16MeV))

approximate_average_flux_16MeV = np.average(all_fluxes_16MeV, weights=1/np.square(all_unc_fluxes_16MeV))
unc_approximate_average_flux_16MeV = np.average(all_unc_fluxes_16MeV)
# print('approximate_average_flux_16MeV: ',approximate_average_flux_16MeV)


# print(y_number_of_atoms)
# print(in_number_of_atoms)
# print(al_number_of_atoms)
# print(zr_number_of_atoms)

production_rate = np.hstack((y_production_rate, in_production_rate, al_production_rate, zr_production_rate))
number_of_atoms = np.hstack((y_number_of_atoms, in_number_of_atoms, al_number_of_atoms, zr_number_of_atoms))
flux_avg_cross_section = np.hstack((y_flux_avg_cross_section, in_flux_avg_cross_section, al_flux_avg_cross_section, zr_flux_avg_cross_section))
unc_production_rate = np.hstack((y_unc_production_rate, in_unc_production_rate, al_unc_production_rate, zr_unc_production_rate))
unc_number_of_atoms = np.hstack((y_unc_number_of_atoms, in_unc_number_of_atoms, al_unc_number_of_atoms, zr_unc_number_of_atoms))
unc_flux_avg_cross_section = np.hstack((y_unc_flux_avg_cross_section, in_unc_flux_avg_cross_section, al_unc_flux_avg_cross_section, zr_unc_flux_avg_cross_section))

print('production_rate',production_rate)
print('number_of_atoms',number_of_atoms)
print('flux_avg_cross_section',flux_avg_cross_section)
print('unc_production_rate',unc_production_rate)
print('unc_number_of_atoms',unc_number_of_atoms)
print('unc_flux_avg_cross_section',unc_flux_avg_cross_section)

# print(len(number_of_atoms))
# plt.errorbar(x_index_33MeV, all_fluxes_33MeV, yerr=all_unc_fluxes_33MeV, marker='.', linestyle='')
# plt.plot(x_index_33MeV, approximate_average_flux_33MeV*np.ones(len(x_index_33MeV)), color='red')
# plt.show()
# plt.errorbar(x_index_16MeV, all_fluxes_16MeV, yerr=all_unc_fluxes_16MeV, marker='.', linestyle='')
# plt.plot(x_index_16MeV, approximate_average_flux_16MeV*np.ones(len(x_index_16MeV)), color='red')
# plt.show()

output_flux, output_unc_flux = Average_Neutron_Flux(production_rate, number_of_atoms, flux_avg_cross_section, unc_production_rate, unc_number_of_atoms, unc_flux_avg_cross_section, csv_filename='averaged_fluxes.csv')

print('all_fluxes_16MeV',all_fluxes_16MeV)
print('all_fluxes_33MeV',all_fluxes_33MeV)


print('approximate_average_flux_16MeV: ',approximate_average_flux_16MeV)
print('approximate_average_flux_33MeV: ',approximate_average_flux_33MeV)

print('correct_average_flux_16MeV: ',output_flux[0])
print('correct_average_flux_33MeV: ',output_flux[1])


    # return output_flux[::-1], output_unc_flux[::-1 ] #returning reversed lists
