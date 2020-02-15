from npat import DecayChain
import csv
import numpy as np

#t_irradiation = 1.0     # in hours
name_of_csv_file = './112In_in33MeV_349.dat'
name_of_csv_file_parent = './112mIn_in33MeV_349.dat'

### 66Ni/66Cu decay chain, units of hours, assumed 2:1 production rate, for t_irradiation hours
#
dc =      DecayChain('112INm', 'h', R={'112INg':1.0, '112INm':10.0}, time=530.0/3600.0)
dc.append(DecayChain('112INm', 'h',                             time=((5.0*60.0)+40.0)/3600.0))
dc.append(DecayChain('112INm', 'h', R={'112IN':1.0, '112INm':10.0}, time=6750.0/3600.0))


### 58m/gCo decay chain, units of hours, assumed 2:1 production rate, for t_irradiation hours
#dc = DecayChain('58COm', 'h', R={'58CO':1.0, '58COm':2.0}, time=t_irradiation)
#dc = DecayChain('66NI', 'h', R={'58CO':1.0, '58COm':2.0}, time=t_irradiation)

# dc.plot()


### Read in numbers of decays from csv file
def decomment(csvfile):
    for row in csvfile:
        raw = row.split('#')[0].strip()
        if raw: yield raw

results = []
with open(name_of_csv_file) as csvfile:
    reader = csv.reader(decomment(csvfile))
    for row in reader:
        results.append(row)

results = np.asarray(results, dtype=float)

results_parent = []
with open(name_of_csv_file_parent) as csvfile:
    reader = csv.reader(decomment(csvfile))
    for row in reader:
        results_parent.append(row)

results_parent = np.asarray(results_parent, dtype=float)



# Reshape csv data structure
list_of_counts = results[:, np.r_[0, 0, 3, 4]]
list_of_counts[:,1] = list_of_counts[:,1] + (results[:,5] / 3600)
print(list_of_counts)

list_of_counts_parent = results_parent[:, np.r_[0, 0, 3, 4]]
list_of_counts_parent[:,1] = list_of_counts_parent[:,1] + (results_parent[:,5] / 3600)
print(list_of_counts_parent)

### Calculate decay over timespan of all counts
dc.append(DecayChain('112INm', 'h', time=max(list_of_counts[:,1])))
# dc.plot()

### Measured counts: [start_time (d), stop_time (d), decays, unc_decays]
### Times relative to last appended DecayChain, i.e. EoB time
dc.counts = {'112INg':list_of_counts, '112INm':list_of_counts_parent}



### Find the scaled production rate that gives us these counts
dc.fit_R( unc=True)
### Only plot the 5 most active isotopes in the decay chain
dc.plot(N_plot=5)
print('              ', dc.isotopes)

EoB_activities = dc.A0

print('Activity (Bq):',EoB_activities)
dc.fit_A0( unc=True)
uncertainty_EoB_activities = dc._unc_A0_fit
print('Uncertainty in Activity (Bq):', uncertainty_EoB_activities)

#print(type(EoB_activities))
#print(type(uncertainty_EoB_activities))

outfile_parent = np.array((EoB_activities[0], uncertainty_EoB_activities[0]))
outfile_daughter = np.array((EoB_activities[1], uncertainty_EoB_activities[1]))

#outfile = np.transpose(np.column_stack((EoB_activities, uncertainty_EoB_activities)))

#print(outfile_parent)
np.savetxt("../find activity/activity_csv/In_112mIn.csv", np.transpose(outfile_parent), delimiter=",")
np.savetxt("../find activity/activity_csv/In_112In.csv", np.transpose(outfile_daughter), delimiter=",")
