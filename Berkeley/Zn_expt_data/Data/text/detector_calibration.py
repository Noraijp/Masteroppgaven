import npat
# from npat import Spectrum, Calibration
import scipy.io
import numpy as np


cb_sp = npat.Spectrum('./AB060818_Ba133_10cm.Spe')
# cb_sp.plot()

cb_sp.meta = {'istp':['133BA'],#, '133BA', '137CS', '60CO'],
				'A0':[39.89E3],#, 38.59E3, 36.7E3, 38.52E3],
				'ref_date':'01/01/2009 12:00:00'}

cb_sp2 = npat.Spectrum('./AC060818_Eu152_10cm.Spe')
# cb_sp.plot()

cb_sp2.meta = {'istp':['152EU'],#, '133BA', '137CS', '60CO'],
				'A0':[39.29E3],#, 38.59E3, 36.7E3, 38.52E3],
				'ref_date':'01/01/2009 12:00:00'}

cb_sp3 = npat.Spectrum('./AA060818_Cs137_10cm.Spe')
# cb_sp.plot()

cb_sp3.meta = {'istp':['137CS'],#, '133BA', '137CS', '60CO'],
				'A0':[38.55E3],#, 38.59E3, 36.7E3, 38.52E3],
				'ref_date':'01/01/2009 12:00:00'}


cb_sp.fit_config = {'xrays':True, 'E_min':45.0, 'bg_fit':True}
# cb_sp.auto_calibrate()
# cb_sp.summarize()

### Plot ADC channels instead of energy
# cb_sp.plot(xcalib=False)

### Pick out a few peaks for manual calibration
# cb_data = [[1890, 356.0]]#,
                        # [1338.5, 244.7],
                        # [1882.5, 344.3],
                        # [2428, 444],
                        # [7698, 1408]]
cb_sp3.fit_config = {'xrays':True, 'E_min':45.0, 'bg_fit':True}
cb_sp2.fit_config = {'xrays':True, 'E_min':45.0, 'bg_fit':True, 'I_min':0.1}

# cb_sp.auto_calibrate(data=cb_data)


# cb_sp2.plot()

cb = npat.Calibration()
# cb_sp.cb = cb
# cb_sp.cb.plot()


sp = npat.Spectrum('./AB060818_Ba133_10cm.Spe')
#sp = npat.Spectrum('/media/andrew/Storage/Documents/School Work/Isotope Data/Ekeberg_Ir_2019/Masteroppgaven/Calibration/room131/text_Spe/BE20190206_Ba133_30cm_room131.Spe')
# sp.meta = {'istp':['58CO','40K','210BI','212BI','210PB','212PB','226RA']}
sp.meta = {'istp':['133BA']}

prev_calib = list(sp.cb.engcal)
sp.fit_config = {'xrays':True,'E_min':45.0, 'bg_fit':True}
# , 'bg_fit':True}
#
# cb = npat.Calibration()
cb.calibrate([cb_sp, cb_sp3, cb_sp2,sp])
sp.auto_calibrate()
exp_engcal = list(sp.cb.engcal)
cb.plot()
#
sp.summarize()
sp.plot()




print(sp.cb.eff(120.0))
print(sp.cb.unc_eff(120.0))

# print("Functional efficiency parameters:")
# print(sp.cb.effcal)

# print("Functional efficiency covariance matrix:")
# print(sp.cb.unc_effcal)


scipy.io.savemat('eff_room131_10.mat', {'popt_10': sp.cb.effcal,'pcov_10': sp.cb.unc_effcal})
