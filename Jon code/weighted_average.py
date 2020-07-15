import numpy as np
import csv
import matplotlib.pyplot as plt
from scipy.constants import elementary_charge



#from des19_BeamCurrent import BeamCurrent




# Nora:
number_of_monitor_foils = 3
monitor_reactions_per_foil = np.array([1, 2, 2])
# Me:
# number_of_monitor_foils = 2
#monitor_reactions_per_foil = np.array([2, 2])
#number_of_monitor_foils = len(monitor_reactions_per_foil)

### Read in numbers of decays from csv file

def Average_Neutron_Flux(production_rate, number_of_atoms, flux_avg_cross_section, unc_production_rate, unc_number_of_atoms, unc_flux_avg_cross_section, solid_angle, csv_filename='averaged_currents.csv'):


    def decomment(csvfile):
	    for row in csvfile:
	        raw = row.split('#')[0].strip()
	        if raw: yield raw

    def read_csv(name_of_csv_file):
        results = []
        with open(name_of_csv_file) as csvfile:
            reader = csv.reader(decomment(csvfile))
            for row in reader:
                results.append(row)

        return np.asarray(results, dtype=float)


    def neutron_flux(R, N_atoms, flux_avg_xs, solid_angle):
        return (R) / (N_atoms * flux_avg_xs * solid_angle)


	# Numerical partial derivatives
    def dIR(R, N_atoms, flux_avg_xs):
        delta_x = 1E-8 * R
        return ((neutron_flux(R + (delta_x/2), N_atoms, flux_avg_xs, solid_angle) - neutron_flux(R - (delta_x/2), N_atoms, flux_avg_xs, solid_angle )) / delta_x)
    def dIdN_atoms(R, N_atoms, flux_avg_xs):
        delta_x = 1E-8 * N_atoms
        return ((neutron_flux(R, N_atoms + (delta_x/2), flux_avg_xs, solid_angle) - neutron_flux(R, N_atoms - (delta_x/2), flux_avg_xs, solid_angle)) / delta_x)
    def dIdXS(R, N_atoms, flux_avg_xs):
        delta_x = 1E-8 * flux_avg_xs
        return ((neutron_flux(R, N_atoms, flux_avg_xs + (delta_x/2), solid_angle) - neutron_flux(R, N_atoms, flux_avg_xs - (delta_x/2), solid_angle)) / delta_x)


	# Approximate uncertainties in neutron flux
    def sigma_flux_approximate(R, N_atoms, flux_avg_xs, unc_R, unc_N_atoms, unc_flux_avg_xs):
        approx_error = np.sqrt(np.power(dIR(R, N_atoms, flux_avg_xs) * unc_R,2) +
        np.power(dIdN_atoms(R, N_atoms, flux_avg_xs) * unc_N_atoms,2) +
        np.power(dIdXS(R, N_atoms, flux_avg_xs) * unc_flux_avg_xs,2))
		# approx_error = 0
		# approx_error = np.power(dIdA0(A0, rho_dr, lambdas, t_irradiation, reaction_integral),2)
        return approx_error


    number_of_monitor_reactions = 0
    #print('yo')
    submatrix_lower_indices = np.zeros(number_of_monitor_foils)
    submatrix_upper_indices = np.zeros(number_of_monitor_foils)


	#### PARAMETERS IN THE FUNCTION des19_BeamCurrent.py
	#A0, dA0, mass_density, sigma_mass_density, lambda_, reaction_integral, uncertainty_integral, irr_time, sigma_irr_time



	# Load in activation data
	#
	# All in nuclei / cm^2
	#                       Cu       Ti
    N_atoms = number_of_atoms
    unc_N_atoms = unc_number_of_atoms

	# All in Bq
	#                      Sc46    V48    Zn62    Zn63
    R = production_rate
    unc_R = unc_production_rate
	# Normalized integral(sigma * dPhidE * dE)
	#
    flux_avg_xs = flux_avg_cross_section
    #unc_rxn_integral = uncertainty_integral #rxn_int * 	percent_rn_uncertainties      #rxn = reactions
    unc_flux_avg_xs = unc_flux_avg_cross_section#rxn_int * 	percent_rn_uncertainties      #rxn = reactions



    for i in range(0, number_of_monitor_foils):
        submatrix_lower_indices[i] = number_of_monitor_reactions
        number_of_monitor_reactions += monitor_reactions_per_foil[i]
        submatrix_upper_indices[i] = number_of_monitor_reactions

    submatrix_lower_indices=submatrix_lower_indices.astype(int)
    submatrix_upper_indices=submatrix_upper_indices.astype(int)


	# Set up correlation matrices
	# Areal density is completely uncorrelated, except within one foil's submatrix
    corr_N_atoms = np.zeros((number_of_monitor_reactions,number_of_monitor_reactions))
	# flux_avg_xs is completely correlated (same for all foils)
    corr_flux_avg_xs = 0.3 * np.ones((number_of_monitor_reactions,number_of_monitor_reactions))
	# production rates are partially uncorrelated (similar subset of efficiencies)
    corr_R = 0.3 * np.ones((number_of_monitor_reactions,number_of_monitor_reactions))    #just set to 0.3 since we do not have MC simulations


	# Set up lists to hold output data
    output_foil_index = []
    output_flux = []
    output_unc_flux = []
    output_percent_unc = []


	# Get correlation submtarix for each monitor foil - n x n, where n= # of reactions per foil
    for i in range(0, number_of_monitor_foils): # Monitor reactions per foil [3,3,1]
        submatrix = np.ones((monitor_reactions_per_foil[i], monitor_reactions_per_foil[i]))
        corr_N_atoms[submatrix_lower_indices[i]:submatrix_upper_indices[i], submatrix_lower_indices[i]:submatrix_upper_indices[i]] = submatrix
        # corr_reaction_integral[submatrix_lower_indices[i]:submatrix_upper_indices[i], submatrix_lower_indices[i]:submatrix_upper_indices[i]] = 0.3*submatrix


	# Ensure all diagonal elements are still ones
    np.fill_diagonal(corr_N_atoms,1)
    np.fill_diagonal(corr_R,1)
    np.fill_diagonal(corr_flux_avg_xs,1)

	# print("corr_lambda")
	# print(corr_lambda)
	# print("corr_areal_density")
	# print(corr_areal_density)
	# print("corr_reaction_integral")
	# print(corr_reaction_integral)
	# print("corr_t_irradiation")
	# print(corr_t_irradiation)
	# print("corr_EoB_activities")
	# print(corr_EoB_activities)


	# Loop over all beam positions
    number_of_energies = len(N_atoms)
    print('number_of_energies: ',number_of_energies)
	# Test mode!!!!
	# number_of_energies = 1

	# Hold curents as we go along...
    fluxes = np.zeros((number_of_energies,number_of_monitor_reactions))
    unc_fluxes = np.zeros((number_of_energies,number_of_monitor_reactions))
	# function_dictionary = {'dIdA0':dIdA0, 'dIdRhoDr':dIdRhoDr, 'dIdLambda':dIdLambda, 'dIdTIrradiation':dIdTIrradiation, 'dIdIntegral':dIdIntegral}
    function_dictionary = {'0':dIR, '1':dIdN_atoms, '2':dIdXS}


    # print('ad: ',areal_density)
    # print('unc_ad: ',uncertainty_areal_density)
    # print('A0:',EoB_activities)
    # print('unc_A0: ',uncertainty_EoB_activities)
    # print('rxn_int: ',reaction_integral)
    # print('unc_rxn_int: ',unc_rxn_int)
    # print('delta_t: ',t_irradiation)
    # print('unc_delta_t: ',uncertainty_t_irradiation)
    # print('loop_lambdas: ',lambdas)
    # print('uncertainty_lambdas: ',uncertainty_lambdas)

    for i_energy in range(0, number_of_energies):
        print('i_energy: ',i_energy)
		# Get nonzero entries in A0:
        nonzero_indices = np.nonzero(R[i_energy,:])
        loop_N_atoms = N_atoms[i_energy,:]
        loop_unc_N_atoms = unc_N_atoms[i_energy,:]
        loop_R = R[i_energy,:]
        loop_unc_R = unc_R[i_energy,:]
        loop_flux_avg_xs = flux_avg_xs[i_energy,:]
        loop_unc_flux_avg_xs = unc_flux_avg_xs[i_energy,:]
        # # delta_t = t_irradiation[i_energy]
        # unc_delta_t = np.ones(number_of_monitor_reactions) *uncertainty_t_irradiation[i_energy]
        # # unc_delta_t = uncertainty_t_irradiation[i_energy]
		# #percent_rn_uncertainties = np.array([0.051054188386, 0.064100768909, 0.084213661384, 0.04385708308])
		# #unc_rxn_int = rxn_int * 	percent_rn_uncertainties      #rxn = reactions
        # unc_rxn_int = unc_rxn_integral[i_energy,:]
        # loop_lambdas = lambdas
        # uncertainty_lambdas = loop_lambdas * 0.001

        # print('loop_N_atoms: ',loop_N_atoms)
        # print('loop_unc_N_atoms: ',loop_unc_N_atoms)
        # print('loop_R:',loop_R)
        # print('loop_unc_R: ',loop_unc_R)
        # print('loop_flux_avg_xs: ',loop_flux_avg_xs)
        # print('loop_unc_flux_avg_xs: ',loop_unc_flux_avg_xs)
        # print('delta_t: ',delta_t)
        # print('unc_delta_t: ',unc_delta_t)
        # print('loop_lambdas: ',loop_lambdas)
        # print('uncertainty_lambdas: ',uncertainty_lambdas)

        if len(np.transpose(nonzero_indices)) == number_of_monitor_reactions:

			# No nonzero indices!!!
			# Keep normal correlation matrices
            loop_corr_N_atoms = corr_N_atoms
            loop_corr_flux_avg_xs = corr_flux_avg_xs
            loop_corr_R = corr_R
            # loop_corr_t_irradiation = corr_t_irradiation
            # loop_corr_EoB_activities = corr_EoB_activities

        else:
			# Some nonzero indices
			# Find which indices are missing!
            temp3 = np.asarray(nonzero_indices[0])
            temp4 = np.array(range(0, number_of_monitor_reactions))
            disjoint_indices = np.setdiff1d(temp4,temp3,assume_unique=False).tolist()
            # print('disjoint indices: ', disjoint_indices)

			# Delete rows and columns in correlation matries of disjoint indices
			# gen = (x for x in xyz if x not in a)
            if len(disjoint_indices) != 1:
                loop_corr_N_atoms = np.delete(corr_N_atoms,np.array(disjoint_indices),0)
                loop_corr_N_atoms = np.delete(loop_corr_N_atoms,np.array(disjoint_indices),1)
                loop_corr_flux_avg_xs = np.delete(corr_flux_avg_xs,np.array(disjoint_indices),0)
                loop_corr_flux_avg_xs = np.delete(loop_corr_flux_avg_xs,np.array(disjoint_indices),1)
                loop_corr_R = np.delete(corr_R,np.array(disjoint_indices),0)
                loop_corr_R = np.delete(loop_corr_R,np.array(disjoint_indices),1)
                # loop_corr_t_irradiation = np.delete(corr_t_irradiation,np.array(disjoint_indices),0)
                # loop_corr_t_irradiation = np.delete(loop_corr_t_irradiation,np.array(disjoint_indices),1)
                # loop_corr_EoB_activities = np.delete(corr_EoB_activities,np.array(disjoint_indices),0)
                # loop_corr_EoB_activities = np.delete(loop_corr_EoB_activities,np.array(disjoint_indices),1)
            else:
                for disjoint_index in disjoint_indices:
                    loop_corr_N_atoms = np.delete(corr_N_atoms,disjoint_index,0)
                    loop_corr_N_atoms = np.delete(loop_corr_N_atoms,disjoint_index,1)
                    loop_corr_flux_avg_xs = np.delete(corr_flux_avg_xs,disjoint_index,0)
                    loop_corr_flux_avg_xs = np.delete(loop_corr_flux_avg_xs,disjoint_index,1)
                    loop_corr_R = np.delete(corr_R,disjoint_index,0)
                    loop_corr_R = np.delete(loop_corr_R,disjoint_index,1)
                    # loop_corr_t_irradiation = np.delete(corr_t_irradiation,disjoint_index,0)
                    # loop_corr_t_irradiation = np.delete(loop_corr_t_irradiation,disjoint_index,1)
                    # loop_corr_EoB_activities = np.delete(corr_EoB_activities,disjoint_index,0)
                    # loop_corr_EoB_activities = np.delete(loop_corr_EoB_activities,disjoint_index,1)

        # print('neutron_flux inputs: ', A0, ad, loop_lambdas[i_energy,:], delta_t, rxn_int)
        # neutron_flux(R, N_atoms, flux_avg_xs)
        temp_fluxes =  neutron_flux(loop_R, loop_N_atoms, loop_flux_avg_xs, solid_angle)
        fluxes[i_energy, :] =  temp_fluxes
        # print('temp_currents: ', temp_currents)
        # sigma_flux_approximate(R, N_atoms, flux_avg_xs, unc_R, unc_N_atoms, unc_flux_avg_xs):
        unc_temp_fluxes = sigma_flux_approximate(loop_R, loop_N_atoms, loop_flux_avg_xs, loop_unc_R, loop_unc_N_atoms, loop_unc_flux_avg_xs)
        unc_fluxes[i_energy,:] = unc_temp_fluxes

        value_array = np.array([loop_R, loop_N_atoms, loop_flux_avg_xs])
        uncertainty_array = np.array([loop_unc_R, loop_unc_N_atoms, loop_unc_flux_avg_xs])
        correlation_array = np.array([loop_corr_R, loop_corr_N_atoms, loop_corr_flux_avg_xs])


		# Set up covariance matrix for current energy position
        cov = np.zeros((len(nonzero_indices[0]),len(nonzero_indices[0])))

		# NaN handling - replace range(0,number_of_monitor_reactions) with indices of nonzero elements of A0?
		# Fill correlation matrices
        for i_index,i_element in enumerate(nonzero_indices[0]):
		# for i in range(0,number_of_monitor_reactions):
            for j_index,j_element in enumerate(nonzero_indices[0]):
			# for j in range(0, number_of_monitor_reactions):
                for dict_index,dict_key in enumerate(function_dictionary):
                    # print("dict_index: ",dict_index)
                    # print("dict_value: ",function_dictionary[dict_key])
                    # print(type(dict_key))
                    # print(A0[i_element], ad[i_element], loop_lambdas[0,i_element], delta_t[i_element], rxn_int[i_element])
                    dIdxi = function_dictionary[dict_key](loop_R[i_element], loop_N_atoms[i_element], loop_flux_avg_xs[i_element])
                    dIdxj = function_dictionary[dict_key](loop_R[j_element], loop_N_atoms[j_element], loop_flux_avg_xs[j_element])
					# print("dIdx_i: ",dIdxi)
					# print("dIdx_j: ",dIdxj)
					# print("unc_xi: ",uncertainty_array[dict_index,i_element])
					# print("unc_xj: ",uncertainty_array[dict_index,j_element])
					# print("corr_x: ", correlation_array[dict_index,i_index,j_index])
                    cov[i_index,j_index] += dIdxi * uncertainty_array[int(dict_key),i_element] *  correlation_array[int(dict_key),i_index,j_index] *  uncertainty_array[int(dict_key),j_element] * dIdxj

    	# print("Final covariance matrix: \n", cov)
        inverted_covariance = np.linalg.inv(cov)
        numerator = 0.0
        denominator = 0.0

        for i_index,i_element in enumerate(nonzero_indices[0]):
            for j_index,j_element in enumerate(nonzero_indices[0]):
                numerator += temp_fluxes[j_element] * inverted_covariance[i_index,j_index]
                denominator +=  inverted_covariance[i_index,j_index]

        weighted_average_flux = numerator/denominator
        uncertainty_weighted_average_flux = np.sqrt(1.0/denominator)


        print("weighted_average_flux: ",weighted_average_flux, " +/- ",uncertainty_weighted_average_flux, " nA     (", 100*uncertainty_weighted_average_flux/weighted_average_flux ," %)")

    	# Append values for current energy
        output_foil_index.append(i_energy)
        output_flux.append(weighted_average_flux)
        output_unc_flux.append(uncertainty_weighted_average_flux)
        output_percent_unc.append(100*uncertainty_weighted_average_flux/weighted_average_flux)
        print("********************************************************************\n")
        print("Raw currents: \n",fluxes)
        # print("Raw unc_currents: \n",unc_currents)


	#matlab_avg_currents = np.array(  [103.6055,   99.3354,  104.7844,  109.5334,  103.4154,   97.9257,   92.6829,   90.6173,   92.6035,   88.0497,   87.8754,   75.7701,   66.8488,    0.1905])
	#matlab_unc_avg_currents = np.array([3.3891,    2.9333,    3.3763,    3.4887,     3.3309,     3.2648,     3.8734,     2.8113,     2.8976,     3.3152,     3.6999,     2.3559,     2.6791,     0.2861])

	# Save final results to csv
    outfile = np.stack((np.transpose(output_foil_index),np.transpose(output_flux),np.transpose(output_unc_flux),np.transpose(output_percent_unc)), axis=-1)
    #import os
    #path = os.getcwd()
    np.savetxt("./{}".format(csv_filename), outfile, delimiter=",", header="Foil Index, Average Neutron Flux (a.u.), Uncertainty in Average Neutron Flux (a.u.), % Uncertainty")



	# Plot output comparisons
	#
	# plt.gca().set_prop_cycle(None)
    output_foil_index2 = np.array(output_foil_index) - 0.2

	# plt.clf()

    plt.gca()
    plt.errorbar(output_foil_index, output_flux, yerr=output_unc_flux, capsize=10.0, markersize=4.0,  marker='s', ls=' ', color='black', capthick=1.5, linewidth=3.0)
	# plt.errorbar(output_foil_index2, matlab_avg_currents, yerr=matlab_unc_avg_currents, capsize=10.0, capthick=2.0, markersize=8.0,  marker='.', ls=' ',  linewidth=2.0)
    for i in range(len(output_foil_index)):
        plt.errorbar(np.ones(len(fluxes[i,:]))*(0.2+output_foil_index[i]), fluxes[i,:], markersize=4.0,  yerr=unc_fluxes[i,:], capsize=5.0,  marker='.', ls=' ',linewidth=0.5, capthick=0.5, color='red')
    plt.legend(['Average Neutron Flux', 'Individual Neutron Fluxes'],loc='best')
	# plt.legend(['Actual Currents', 'Approximate Currents', 'Individual Currents'],loc='lower left')
    plt.show()


    #output_mu = output_mu.reverse()

    return output_flux, output_unc_flux
    # return output_flux[::-1], output_unc_flux[::-1 ] #returning reversed lists
