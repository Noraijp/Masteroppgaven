import numpy as np
import matplotlib.pyplot as plt
import npat
import csv

from npat import *


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

	# print('\n************************************************\n')
	return flux_array_33MeV, unc_flux_array_33MeV, flux_array_16MeV, unc_flux_array_16MeV






## for yttrium

csv_list      = ['Y_88Y.csv']
reaction_list = ['89Y(n,2n)88Y']
target_list   = ['89Y']
product_list  = ['88Y']
mass_33MeV = 0.4657 #g
mass_16MeV = 0.5053
unc_mass_33MeV = 0.0006  #g
unc_mass_16MeV = 0.0021

y_avg_flux_33MeV, y_unc_avg_flux_33MeV, y_avg_flux_16MeV, y_unc_avg_flux_16MeV = calculate_flux(csv_list, reaction_list, target_list, product_list, mass_16MeV, unc_mass_16MeV, mass_33MeV, unc_mass_33MeV)
# print('Y flux at 33MeV           : ', y_avg_flux_33MeV)
# print('Y flux uncertinty_33MeV : ', y_unc_avg_flux_33MeV)
# print('Y flux at 16MeV           : ', y_avg_flux_16MeV)
# print('Y flux uncertinty_16MeV : ', y_unc_avg_flux_16MeV)



### for indium

csv_list      = ['In_114mIn.csv','In_114mIn.csv' ,'In_114mIn.csv', 'In_113mIn.csv',  'In_115mIn.csv',  'In_116mIn.csv']
#reaction_list = ['113IN(n,g)114INm', '113IN(n,inl)113INm', '115IN(n,inl)115INm', '115IN(n,g)116INm']
reaction_list = ['natIN(n,x)114INm', '113IN(n,g)114INm', '115IN(n,2n)114INm', '113IN(n,inl)113INm', '115IN(n,inl)115INm', '115IN(n,g)116INm']
target_list   = ['natIN', '113IN', '115IN', '113IN', '115IN', '115IN']
product_list  = ['114INm', '114INm', '114INm','113INm', '115INm', '116INm']
mass_33MeV = 0.5443 #g
unc_mass_33MeV = 0.0013  #g
mass_16MeV = 0.5530 #g
unc_mass_16MeV = 0.0035  #g

# ## Search the IRDFF neutron library for reactions
# lb = Library('IRDFF')
# print(lb.search(target='113IN'))
# print(lb.search(product='114INm'))
# #print(lb.search(target='113IN',product='114INm1'))

in_avg_flux_33MeV, in_unc_avg_flux_33MeV, in_avg_flux_16MeV, in_unc_avg_flux_16MeV = calculate_flux(csv_list, reaction_list, target_list, product_list, mass_16MeV, unc_mass_16MeV, mass_33MeV, unc_mass_33MeV)
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

al_avg_flux_33MeV, al_unc_avg_flux_33MeV, al_avg_flux_16MeV, al_unc_avg_flux_16MeV = calculate_flux(csv_list, reaction_list, target_list, product_list, mass_16MeV, unc_mass_16MeV, mass_33MeV, unc_mass_33MeV)
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

zr_avg_flux_33MeV, zr_unc_avg_flux_33MeV, zr_avg_flux_16MeV, zr_unc_avg_flux_16MeV = calculate_flux(csv_list, reaction_list, target_list, product_list, mass_16MeV, unc_mass_16MeV, mass_33MeV, unc_mass_33MeV)
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
print(all_fluxes_33MeV)
# print(all_fluxes[3])
x_index_33MeV = range(len(all_fluxes_33MeV))

approximate_average_flux_33MeV = np.average(all_fluxes_33MeV, weights=1/np.square(all_unc_fluxes_33MeV))
unc_approximate_average_flux_33MeV = np.average(all_unc_fluxes_33MeV)
print('approximate_average_flux_33MeV: ',approximate_average_flux_33MeV)


all_fluxes_16MeV = []
all_unc_fluxes_16MeV = []

for i in (y_avg_flux_16MeV, in_avg_flux_16MeV, al_avg_flux_16MeV, zr_avg_flux_16MeV):
	all_fluxes_16MeV.extend(i[:])
for i in (y_unc_avg_flux_16MeV, in_unc_avg_flux_16MeV, al_unc_avg_flux_16MeV, zr_unc_avg_flux_16MeV):
	all_unc_fluxes_16MeV.extend(i[:])
print(all_fluxes_16MeV)
# print(all_fluxes[3])
x_index_16MeV = range(len(all_fluxes_16MeV))

approximate_average_flux_16MeV = np.average(all_fluxes_16MeV, weights=1/np.square(all_unc_fluxes_16MeV))
unc_approximate_average_flux_16MeV = np.average(all_unc_fluxes_16MeV)
print('approximate_average_flux_16MeV: ',approximate_average_flux_16MeV)

# plt.errorbar(x_index_33MeV, all_fluxes_33MeV, yerr=all_unc_fluxes_33MeV, marker='.', linestyle='')
# plt.plot(x_index_33MeV, approximate_average_flux_33MeV*np.ones(len(x_index_33MeV)), color='red')
# plt.show()
# plt.errorbar(x_index_16MeV, all_fluxes_16MeV, yerr=all_unc_fluxes_16MeV, marker='.', linestyle='')
# plt.plot(x_index_16MeV, approximate_average_flux_16MeV*np.ones(len(x_index_16MeV)), color='red')
# plt.show()


# def calculate_flux(csv_list, reaction_list, target_list, product_list, mass_16MeV, unc_mass_16MeV, mass_33MeV, unc_mass_33MeV):
def calculate_xs(product_name, target_name, csv_list, mass_16MeV, mass_33MeV, unc_mass_16MeV, unc_mass_33MeV):

	### for zinc
	# product_name='67CU'
	# itp = npat.Isotope(target_name[j])
	xs_array_16MeV = []
	xs_array_33MeV = []


	# product_activity_data = read_csv('Zn_67Cu.csv')
	#print('CSV data: \n',product_activity_data)
	for j in range(len(product_name)):
		itp = npat.Isotope(target_name[j])
		product_activity_data = read_csv(csv_list[j])
		product_activity_16MeV = product_activity_data[:,0]
		product_activity_33MeV = product_activity_data[:,1]
		# print('product_activity_16MeV', product_activity_16MeV[0])

		product_dc_33MeV =      DecayChain(product_name[j], R=1, time=530.0)
		product_dc_33MeV.append(DecayChain(product_name[j],      time=((5.0*60.0)+40.0)))
		product_dc_33MeV.append(DecayChain(product_name[j], R=1, time=6750.0))

		product_dc_16MeV =      DecayChain(product_name[j], R=1, time= (7*3600) + (57*60))

		#dc = npat.DecayChain(product_name, R=1, time=irr_time)
		product_dc_33MeV.append(npat.DecayChain(product_name[j], time=1E6))
		product_dc_16MeV.append(npat.DecayChain(product_name[j], time=1E6))

		count_start_time = 0.0 # h
		#measured_A0_67ZN_activity = 3.7E3 # bq

		R_33MeV = product_activity_33MeV[0]/product_dc_33MeV.activity(product_name[j], t=count_start_time, units='h')
		R_16MeV = product_activity_16MeV[0]/product_dc_16MeV.activity(product_name[j], t=count_start_time, units='h')
		# R_33MeV = product_activity_33MeV[0]
		# R_16MeV = product_activity_16MeV[0]


		# itp = npat.Isotope('natZN')
		if str(itp)[0:3] == 'nat':
			ab = 1.0
		else:
			ab = 1E-2*itp.abundance()

		if str(itp) == 'natIN':
			isotope_mass = 114.818
		elif str(itp) == 'natY':
			isotope_mass = 88.90584
		elif str(itp) == 'natAL':
			isotope_mass = 26.9815384
		elif str(itp) == 'natZR':
			isotope_mass = 91.224
		elif str(itp) == 'natZN':
			isotope_mass = 65.38
		else:
			isotope_mass = itp.mass

		n_atoms_33MeV = (ab*mass_33MeV*6.022E-1)/isotope_mass
		n_atoms_16MeV = (ab*mass_16MeV*6.022E-1)/isotope_mass


		xs_33MeV = R_33MeV/(n_atoms_33MeV*approximate_average_flux_33MeV)
		xs_16MeV = R_16MeV/(n_atoms_16MeV*approximate_average_flux_16MeV)
		print('\n************************************************\n')
		print('For reaction ',target_name[j],'(n,x)',product_name[j],':')
		print('Our measured 33MeV cross section (mb); ', xs_33MeV)
		print('Our measured 16MeV cross section (mb); ', xs_16MeV)

		xs_array_33MeV.append(xs_33MeV)
		xs_array_16MeV.append(xs_16MeV)

		lb = Library('IRDFF')
		#print(lb.search(target='113IN'))
		IRDFF_reactions_found = lb.search(product=product_name[j])
		for rxn in IRDFF_reactions_found:
			# print(rxn)
			length_target_string  = len(target_name[j])
			length_product_string = len(product_name[j])
			# print(rxn[:length_target_string])
			if rxn[:length_target_string] == target_name[j]:
				if rxn[-length_product_string:] == product_name[j]:
					# print(rxn[:length_target_string])
					# print(rxn[-length_product_string:])
					# we are looking at the correct reaction channel for IRDFF!

					rx = npat.Reaction(rxn, library='IRDFF')
					meulders_33MeV_xs_avg = rx.average(meulders_33MeV[:,0], meulders_33MeV[:,1])
					meulders_16MeV_xs_avg = rx.average(meulders_16MeV[:,0], meulders_16MeV[:,1])
					print('meulders_33MeV_xs_avg: ',meulders_33MeV_xs_avg)
					print('meulders_16MeV_xs_avg: ',meulders_16MeV_xs_avg)

	return xs_array_16MeV, xs_array_33MeV



### zink

product_name  = ['67CU', '67CU', '64CU', '64CU', '62ZN', '63ZN', '65NI', '65ZN', '66CU', '66NI', '69ZNm']
target_name  = ['natZN','67ZN', 'natZN', '64ZN', 'natZN', 'natZN', 'natZN', 'natZN', 'natZN', 'natZN', 'natZN']
csv_list = ['Zn_67Cu.csv', 'Zn_67Cu.csv', 'Zn_64Cu.csv', 'Zn_64Cu.csv', 'Zn_62Zn.csv', 'Zn_63Zn.csv',
            'Zn_65Ni.csv','Zn_65Zn.csv', 'Zn_66Cu.csv', 'Zn_66Ni.csv', 'Zn_69mZn.csv']
mass_33MeV = 0.8427 #g
unc_mass_33MeV = 0.0006
mass_16MeV = 0.8463 #g
unc_mass_16MeV = 0.0015

zn_xs_16MeV, zn_xs_33MeV = calculate_xs(product_name, target_name, csv_list, mass_16MeV, mass_33MeV, unc_mass_16MeV, unc_mass_33MeV)




### Zirconium

product_name = ['90Ym', '91Ym', '91SR', '92Y', '93Y', '95NB', '95ZR', '97NB', '97NB', '97ZR', '98ZR']
target_name = ['natZR', 'natZR', 'natZR', 'natZR', 'natZR', 'natZR', 'natZR', 'natZR', 'natZR', 'natZR', 'natZR']
csv_list = ['Zr_90mY.csv', 'Zr_91mY.csv', 'Zr_91Sr.csv', 'Zr_92Y.csv', 'Zr_93Y.csv', 'Zr_95Nb.csv',
            'Zr_95Zr.csv', 'Zr_97Nb_33.csv', 'Zr_97Nb.csv', 'Zr_97Zr.csv', 'Zr_98Zr.csv']
mass_33MeV = 0.7557 #g
unc_mass_33MeV = 0.0012  #g
mass_16MeV = 0.7560 #g
unc_mass_16MeV = 0.0010 #g

zn_xs_16MeV, zn_xs_33MeV = calculate_xs(product_name, target_name, csv_list, mass_16MeV, mass_33MeV, unc_mass_16MeV, unc_mass_33MeV)



# ### for indium

product_name = ['111IN', '112IN', '112INm', ]
target_lname = ['natIN', 'natIN', 'natIN']
csv_list = ['In_111In.csv', 'In_112In.csv', 'In_112mIn.csv']
mass_33MeV = 0.5443 #g
unc_mass_33MeV = 0.0013  #g
mass_16MeV = 0.5530 #g
unc_mass_16MeV = 0.0035  #g

zn_xs_16MeV, zn_xs_33MeV = calculate_xs(product_name, target_name, csv_list, mass_16MeV, mass_33MeV, unc_mass_16MeV, unc_mass_33MeV)



# ## for yttrium
csv_list      = ['Y_88Y.csv']
product_name = ['87SRm', '87Ym', '87Y', '90Ym']
target_name = ['natY', 'natY', 'natY', 'natY']
csv_list = ['Y_87mSr.csv', 'Y_87mY.csv', 'Y_87Y.csv', 'Y_90mY.csv']
mass_33MeV = 0.4657 #g
mass_16MeV = 0.5053
unc_mass_33MeV = 0.0006  #g
unc_mass_16MeV = 0.0021

zn_xs_16MeV, zn_xs_33MeV = calculate_xs(product_name, target_name, csv_list, mass_16MeV, mass_33MeV, unc_mass_16MeV, unc_mass_33MeV)



### Aluminum
#
# product_name = ['24NA', '24NA']
# target_name = ['27AL', '27AL']
# csv_list = ['Al_24Na.csv', 'Al_24Na.csv']
# mass_33MeV = 0.2563 #g
# unc_mass_33MeV = 0.0015  #g
# mass_16MeV = 0.2573 #g
# unc_mass_16MeV = 0.0006 #g
#
# zn_xs_16MeV, zn_xs_33MeV = calculate_xs(product_name, target_name, csv_list, mass_16MeV, mass_33MeV, unc_mass_16MeV, unc_mass_33MeV)
