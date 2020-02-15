import numpy as np
import matplotlib.pyplot as plt
import npat

def calc_SA(R_beam=0.5, R_sample=0.05, dist=5.0, N=1E6):
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

print(calc_SA(2.5, 5.6, 5.0)/(4.0*np.pi))

R_sample = [0.05, 0.5, 1.0, 2.0, 5.0, 10.0, 100.0]
SA = []
for r in R_sample:
	SA.append(calc_SA(R_sample=r))

plt.plot(R_sample, SA)
plt.show()