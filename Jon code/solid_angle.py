import numpy as np
import matplotlib.pyplot as plt
import npat

def calc_SA(R_beam=0.5, R_sample=0.05, dist=5.0, N=1E8):
	#### Solid angle between two concentric circles
	N = int(N)
	r = np.sqrt(np.random.uniform(size=N))*R_beam
	ang = np.pi*2.0*np.random.uniform(size=N)

	x0 = r*np.cos(ang)
	y0 = r*np.sin(ang)

	tht = np.arccos(2.0*np.random.uniform(size=N)-1.0)
	phi = 2.0*np.pi*np.random.uniform(size=N)

	xf = x0 + dist*np.tan(tht)*np.cos(phi)
	yf = y0 + dist*np.tan(tht)*np.sin(phi)

	return 2.0*np.pi*float(len(np.where(np.sqrt(xf**2+yf**2)<=R_sample)[0]))/float(N)

# print(calc_SA(2.5, 5.6, 5.0)/(4.0*np.pi))

# R_sample = [0.05, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0]
# SA = []
# for r in R_sample:
# 	SA.append(calc_SA(R_sample=r))

# plt.plot(R_sample, SA)
# plt.show()

#REDIGER det skal være radiusen til beamen som ble kjørt før eksperimentet
#the radius (in cm) of the deuteron beam, which we know teh size of from your gafchromic films
R_beam_16MeV = 0.663324958 #apporox
R_beam_33MeV = 1.100843454 #approx

output_foil_index = np.array([0,1])





R_sample_zn_16MeV = 0.637 #cm
R_sample_zn_33MeV = 0.644
dist_to_zn_16MeV = 10.1 #cm
dist_to_zn_33MeV = 10.1

#zn_solid_angle = np.array([calc_SA(R_beam=R_beam_16MeV, R_sample=R_sample_zn_16MeV, dist=dist_to_zn_16MeV), calc_SA(R_beam=R_beam_33MeV, R_sample=R_sample_zn_33MeV, dist=dist_to_zn_33MeV)])
#print('Zn Solid Angle: ',zn_solid_angle)
#outfile = np.stack((np.transpose(output_foil_index),np.transpose(zn_solid_angle)), axis=-1)
#np.savetxt("./{}".format('zn_solid_angle.csv'), outfile, delimiter=",", header="Foil Index, Solid Angle")



R_sample_y_16MeV = 0.633 #cm
R_sample_y_33MeV = 0.590
dist_to_y_16MeV = 10.7 #cm
dist_to_y_33MeV = 11.1

#y_solid_angle = np.array([calc_SA(R_beam=R_beam_16MeV, R_sample=R_sample_y_16MeV, dist=dist_to_y_16MeV), calc_SA(R_beam=R_beam_33MeV, R_sample=R_sample_y_33MeV, dist=dist_to_y_33MeV)])
#print('Y Solid Angle: ',y_solid_angle)
#outfile = np.stack((np.transpose(output_foil_index),np.transpose(y_solid_angle)), axis=-1)
#np.savetxt("./{}".format('y_solid_angle.csv'), outfile, delimiter=",", header="Foil Index, Solid Angle")




R_sample_al_16MeV = 0.6442 #cm
R_sample_al_33MeV = 0.646
dist_to_al_16MeV = 10.5 #cm
dist_to_al_33MeV = 10.5

#al_solid_angle = np.array([calc_SA(R_beam=R_beam_16MeV, R_sample=R_sample_al_16MeV, dist=dist_to_al_16MeV), calc_SA(R_beam=R_beam_33MeV, R_sample=R_sample_al_33MeV, dist=dist_to_al_33MeV)])
#print('Al Solid Angle: ',al_solid_angle)
#outfile = np.stack((np.transpose(output_foil_index),np.transpose(al_solid_angle)), axis=-1)
#np.savetxt("./{}".format('al_solid_angle.csv'), outfile, delimiter=",", header="Foil Index, Solid Angle")



R_sample_zr_16MeV = 0.615 #cm
R_sample_zr_33MeV = 0.615
dist_to_zr_16MeV = 10.9 #cm
dist_to_zr_33MeV = 10.8

#zr_solid_angle = np.array([calc_SA(R_beam=R_beam_16MeV, R_sample=R_sample_zr_16MeV, dist=dist_to_zr_16MeV), calc_SA(R_beam=R_beam_33MeV, R_sample=R_sample_zr_33MeV, dist=dist_to_zr_33MeV)])
#print('Zr Solid Angle: ',zr_solid_angle)
#outfile = np.stack((np.transpose(output_foil_index),np.transpose(zr_solid_angle)), axis=-1)
#np.savetxt("./{}".format('zr_solid_angle.csv'), outfile, delimiter=",", header="Foil Index, Solid Angle")


R_sample_in_16MeV = 0.701 #cm
R_sample_in_33MeV = 0.699
dist_to_in_16MeV = 10.3 #cm #### !!!!!!! Doubble check this!! IMPOTANT!!!!!!!
dist_to_in_33MeV = 10.3

in_solid_angle = np.array([calc_SA(R_beam=R_beam_16MeV, R_sample=R_sample_in_16MeV, dist=dist_to_in_16MeV), calc_SA(R_beam=R_beam_33MeV, R_sample=R_sample_in_33MeV, dist=dist_to_in_33MeV)])
print('In Solid Angle: ',in_solid_angle)
outfile = np.stack((np.transpose(output_foil_index),np.transpose(in_solid_angle)), axis=-1)
np.savetxt("./{}".format('in_solid_angle.csv'), outfile, delimiter=",", header="Foil Index, Solid Angle")
