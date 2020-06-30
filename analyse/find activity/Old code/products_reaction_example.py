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


def read_flux_csv(name_of_csv_file):
	results = []
	with open(name_of_csv_file) as csvfile:
	    reader = csv.reader(decomment(csvfile))
	    for row in reader:
	        results.append(row)

	return np.asarray(results, dtype=float)



true_neutron_fluxes =  read_flux_csv('./averaged_fluxes.csv')
neutron_flux_16MeV = true_neutron_fluxes[0,:]
neutron_flux_33MeV = true_neutron_fluxes[1,:]
approximate_average_flux_16MeV = neutron_flux_16MeV[1]
approximate_average_flux_33MeV = neutron_flux_33MeV[1]
unc_approximate_average_flux_16MeV = neutron_flux_16MeV[2]
unc_approximate_average_flux_33MeV = neutron_flux_33MeV[2]


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
		# print('product_activity_33MeV', product_activity_33MeV[0])

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
		# print('% uncertainty in activity; ',100*product_activity_33MeV[1]/product_activity_33MeV[0])
		# print('% uncertainty in flux; ',(100*unc_approximate_average_flux_33MeV/approximate_average_flux_33MeV))
		unc_xs_33MeV = xs_33MeV*np.sqrt(  (product_activity_33MeV[1]/product_activity_33MeV[0])**2  + (unc_approximate_average_flux_33MeV/approximate_average_flux_33MeV)**2 + (unc_mass_33MeV/mass_33MeV)**2   )
		unc_xs_16MeV = xs_16MeV*np.sqrt(  (product_activity_16MeV[1]/product_activity_16MeV[0])**2  + (unc_approximate_average_flux_16MeV/approximate_average_flux_16MeV)**2 + (unc_mass_16MeV/mass_16MeV)**2   )

		# print('\n************************************************\n')
		# print('For reaction ',target_name[j],'(n,x)',product_name[j],':')
		# print('product_activity_16MeV: {0:.2f} Bq' .format(product_activity_16MeV[0]))
		# print('product_activity_33MeV: {0:.2f} Bq' .format(product_activity_33MeV[0]))
		# print('Our measured 16MeV cross section: {0:.2f} +/- {1:.2f} mb ({2:.2f}%)'.format(xs_16MeV,unc_xs_16MeV,100*unc_xs_16MeV/xs_16MeV))
		# print('Our measured 33MeV cross section: {0:.2f} +/- {1:.2f} mb ({2:.2f}%)'.format(xs_33MeV,unc_xs_33MeV,100*unc_xs_33MeV/xs_33MeV))


		E = [16, 33] #MeV
		CS = [xs_16MeV, xs_33MeV]
		dCS = [unc_xs_16MeV, unc_xs_33MeV]

		#save to file
		csv_save_array = np.vstack((E, CS, dCS)).T
		print(csv_save_array)

		path_to_cs_csv = os.getcwd() + '/CrossSections/CrossSections_csv/'
		print(path_to_cs_csv)
        np.savetxt(path_to_cs_csv  + reaction, csv_save_array, delimiter=',', header='E, dE, CS, dCS', fmt="%s"  )#, %.6f, %.6f")

		#gaar inn i moulders
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
					print('meulders_16MeV_xs_avg: {0:.2f} mb'.format(meulders_16MeV_xs_avg))
					print('meulders_33MeV_xs_avg: {0:.2f} mb'.format(meulders_33MeV_xs_avg))

	return xs_array_16MeV, xs_array_33MeV



# # ### zink
#
# product_name  = ['67CU', '64CU', '62ZN', '63ZN', '65NI', '65ZN', '66CU', '66NI', '69ZNm']
# target_name  = ['natZN', 'natZN', 'natZN', 'natZN', 'natZN', 'natZN', 'natZN', 'natZN', 'natZN']
# csv_list = ['Zn_67Cu.csv', 'Zn_64Cu.csv', 'Zn_62Zn.csv', 'Zn_63Zn.csv',
#             'Zn_65Ni.csv','Zn_65Zn.csv', 'Zn_66Cu.csv', 'Zn_66Ni.csv', 'Zn_69mZn.csv']
# mass_33MeV = 0.8427 #g
# unc_mass_33MeV = 0.0006
# mass_16MeV = 0.8463 #g
# unc_mass_16MeV = 0.0015
#
# zn_xs_16MeV, zn_xs_33MeV = calculate_xs(product_name, target_name, csv_list, mass_16MeV, mass_33MeV, unc_mass_16MeV, unc_mass_33MeV)

#
#

### Zirconium

product_name = ['90Ym', '91Ym', '91SR', '92Y', '93Y', '95NB', '95ZR', '97NB_33', '97NB', '97ZR', '98ZR','89ZR']
target_name = ['natZR', 'natZR', 'natZR', 'natZR', 'natZR', 'natZR', 'natZR', 'natZR', 'natZR', 'natZR', 'natZR', 'natZR']
csv_list = ['Zr_90mY.csv', 'Zr_91mY.csv', 'Zr_91Sr.csv', 'Zr_92Y.csv', 'Zr_93Y.csv', 'Zr_95Nb.csv' ,
            'Zr_95Zr.csv', 'Zr_97Nb_33.csv', 'Zr_97Nb.csv', 'Zr_97Zr.csv', 'Zr_98Zr.csv', 'Zr_89Zr.csv']
mass_33MeV = 0.7557 #g
unc_mass_33MeV = 0.0012  #g
mass_16MeV = 0.7560 #g
unc_mass_16MeV = 0.0010 #g

zn_xs_16MeV, zn_xs_33MeV = calculate_xs(product_name, target_name, csv_list, mass_16MeV, mass_33MeV, unc_mass_16MeV, unc_mass_33MeV)


#
# # ### for indium
#
# product_name = ['111IN', '112IN', '112INm']
# target_name = ['natIN', 'natIN', 'natIN']
# csv_list = ['In_111In.csv', 'In_112In.csv', 'In_112mIn.csv']
# mass_33MeV = 0.5443 #g
# unc_mass_33MeV = 0.0013  #g
# mass_16MeV = 0.5530 #g
# unc_mass_16MeV = 0.0035  #g
#
# zn_xs_16MeV, zn_xs_33MeV = calculate_xs(product_name, target_name, csv_list, mass_16MeV, mass_33MeV, unc_mass_16MeV, unc_mass_33MeV)
#
#
#
# # ## for yttrium
# product_name = ['87SRm', '87Ym', '87Y', '90Ym']
# target_name = ['natY', 'natY', 'natY', 'natY']
# csv_list = ['Y_87mSr.csv', 'Y_87mY.csv', 'Y_87Y.csv', 'Y_90mY.csv']
# mass_33MeV = 0.4657 #g
# mass_16MeV = 0.5053
# unc_mass_33MeV = 0.0006  #g
# unc_mass_16MeV = 0.0021
#
# zn_xs_16MeV, zn_xs_33MeV = calculate_xs(product_name, target_name, csv_list, mass_16MeV, mass_33MeV, unc_mass_16MeV, unc_mass_33MeV)
#
#

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
