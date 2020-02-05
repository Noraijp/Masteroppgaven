import numpy as np, matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit

m = 60; h = m*60; d = h*24; y = d*356  #converters to seconds

#fnames = '/Users/Nora/Documents/UiO/Masteroppgaven/analyse/csv/fnames.txt'
#with open (fnames, "r") as myfile:#
#    data=myfile.readlines()
    #print(type(data))


path = '/Users/Nora/Documents/UiO/Masteroppgaven/analyse/csv/'


#comment for testing git

####Zink foils########
def Zn_62Zn(): #check,  #mon, single
    foil1 = path + '62Zn_zn16MeV_130.dat'
    foil2 = path + '62Zn_33MeV_230.dat'
    #foil3 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_329.dat'
    #foil4 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_429.dat'
    #foil5 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_529.dat'
    #foil6 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_629.dat'
    #foil7 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_729.dat'
    #foil8 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_829.dat'
    #foil9 = '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_929.dat'
    #foil10= '/Users/hannah/Documents/UIO/Masteroppgaven/Data/Data_analysis/csv/62Zn_1029.dat'
    list = [foil1, foil2]  #, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_ = np.log(2)/(9.186*h)
    return list, lambda_    #mon, single

def Zn_67Cu(): #check
    foil1 = path + '67Cu_zn16MeV_130.dat'
    foil2 = path + '67Cu_33MeV_230.dat'
    list = [foil1, foil2] #, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_ = np.log(2)/(61.83*h)
    return list, lambda_    #mon, single


def Zn_65Ni(): #peak is wrong
    foil1 = path + '65Ni_zn16MeV_130.dat'
    foil2 = path + '65Ni_zn33MeV_230.dat'
    list = [foil1, foil2] #, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_ = np.log(2)/(2.5175*h) #52mMn
    return list, lambda_


#def Zn_XX(): #two steps  decay
#    foil1 = path + ''
#    foil2 = path + ''
#    list = [foil1, foil2] #, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
#    lambda_parent = np.log(2)/(2.5175*h) #52mMn
#    lambda_daughter = np.log(2)/(halflife_XX * h)
#    return list, lambda_parent, lambda_daughter


def Zn_69mZn():   #single decay        #look wierd
    foil1 = path + '69mZn_zn16MeV_130.dat'
    foil2 = path + '69mZn_zn33MeV_230.dat'
    list = [foil1, foil2]
    #lambda_parent = np.log(2)/(6.075*d)
    lambda_ = np.log(2)/(13.756*h)
    return list, lambda_#parent, lambda_daughter


def Zn_63Zn():   #double decay from 57Ni
    foil1 = path + '63Zn_zn16MeV_130.dat'
    foil2 = path + '63Zn_zn33MeV_230.dat'
    list = [foil1, foil2]
    lambda_ = np.log(2)/(38.47*m)
    return list, lambda_


def Zn_66Cu():  #single decay
    foil1 = path + '66Cu_zn33MeV_230.dat'
    #foil2 = path + ''
    list = [foil1]
    lambda_ = np.log(2)/(5.120*m)
    return list, lambda_


def Zn_64Cu():   #double decay with unkown parent 58mCo   #weird acting
    foil1 = path + '64Cu_zn16MeV_130.dat'
    foil2 = path + '64Cu_zn33MeV_230.dat'
    list = [foil1, foil2]
    lambda_ = np.log(2)/(12.701*h)
    return list, lambda_


def Zn_24Na():   #double decay with unkown parent 58mCo   #weird acting
    foil1 = path + '24Na_zn16MeV_130.dat'
    foil2 = path + '24Na_zn33MeV_230.dat'
    list = [foil1, foil2]
    lambda_ = np.log(2)/(14.997*h)
    return list, lambda_


def Zn_56Mn(): # form iron when we stamped the zink foils
    foil1 = path + '56Mn_Zn16MeV_130.dat'
    list = [foil1]
    lambda_ = np.log(2)/(2.5789*h)
    return list, lambda_


def Zn_40K():
    foil1 = path + '40K_zn16MeV_130.dat'
    list = [foil1]
    lambda_ = np.log(2)/(1.248E9*y)
    return list, lambda_


#####Zirconium FOILS######
#Lacking: 60Co

def Zr_90mY(): #non-mon, single
    foil1 = path + '90mY_zr33MeV_240.dat'
    list = [foil1]
    lambda_ = np.log(2)/(3.19*h)
    return list, lambda_


def Zr_93Y(): #non-mon, single
    foil1 = path + '93Y_zr33MeV_240.dat'
    list = [foil1]
    lambda_ = np.log(2)/(10.18*h)
    return list, lambda_

def Zr_92Y(): #non-mon, single
    foil1 = path + '92Y_zr33MeV_240.dat'
    list = [foil1]
    lambda_ = np.log(2)/(3.54*h)
    return list, lambda_

#def Zr_91mY(): #non-mon, single
#    foil1 = path + '91mY_zr33MeV_240.dat'
#    list = [foil1]
#    lambda_ = np.log(2)/(49.71*m)
#    return list, lambda_


def Zr_91mY(): #two steps decay #Den kjernen du ser på
    foil1 = path + '91mY_zr33MeV_240.dat'
    list = [foil1] #, foil3, foil4, foil5, foil6, foil7, foil8, foil9, foil10]
    lambda_parent = np.log(2)/(9.63*h) #91Sr
    lambda_daughter = np.log(2)/(49.71*m) #91mY
    return list, lambda_parent, lambda_daughter


def Zr_97Nb(): #two step decay
    foil1 = path + '97Nb_zr16MeV_140.dat'
    foil2 = path + '97Nb_zr33MeV_240.dat'
    list = [foil1, foil2]
    lambda_parent = np.log(2)/(16.749*h) #97Zr
    lambda_daughter = np.log(2)/(72.1*m)
    return list, lambda_parent, lambda_daughter


def Zr_95Zr(): #non-mon, single
    foil1 = path + '95Zr_zr16MeV_140.dat'
    foil2 = path + '95Zr_zr33MeV_240.dat'
    list = [foil1, foil2]
    lambda_ = np.log(2)/(64.032*d)
    return list, lambda_


def Zr_97Zr(): #non-mon, single
    foil1 = path + '97Zr_zr16MeV_140.dat'
    list = [foil1]
    lambda_ = np.log(2)/(16.749*h)
    return list, lambda_

def Zr_98Zr(): #non-mon, single
    foil1 = path + '98Zr_zr16MeV_140.dat'
    list = [foil1]
    lambda_ = np.log(2)/(30.7*s)
    return list, lambda_

def Zr_91Sr(): #non-mon, single
    foil1 = path + '91Sr_zr33MeV_240.dat'
    list = [foil1]
    lambda_ = np.log(2)/(9.65*h)
    return list, lambda_

def Zr_89Zr(): #non-mon, single
    foil1 = path + '89Zr_zr33MeV_240.dat'
    list = [foil1]
    lambda_ = np.log(2)/(78.41*h)
    return list, lambda_

def Zr_24Na(): #non-mon, single
    foil1 = path + '24Na_zr33MeV_240.dat'
    list = [foil1]
    lambda_ = np.log(2)/(14.997*h)
    return list, lambda_

def Zr_95Nb(): #non-mon, single
    foil1 = path + '95Nb_zr16MeV_140.dat'
    list = [foil1]
    lambda_daughter = np.log(2)/(34.992*d)
    lambda_parent = np.log(2)/(64.032*d) #Zr95
    return list, lambda_parent, lambda_daughter



    ### YTTRIUM FOILS ###

def Y_90mY(): #non-mon, single
    foil1 = path + '90mY_Y33MeV_239.dat'
    list = [foil1]
    lambda_ = np.log(2)/(3.19*h)
    return list, lambda_


def Y_87mY(): #non-mon, single
    foil1 = path + '87mY_Y33MeV_239.dat'
    list = [foil1]
    lambda_ = np.log(2)/(13.37*h)
    return list, lambda_


def Y_87Y(): #double decay
    foil1 = path + '87Y_Y33MeV_239.dat'
    list = [foil1]
    lambda_parent = np.log(2)/(3.19*h) # 87mY
    lambda_daughter = np.log(2)/(79.8*h)
    return list, lambda_parent, lambda_daughter

def Y_88Y(): #non-mon, single
    foil1 = path + '88Y_Y16MeV_139.dat'
    foil2 = path + '88Y_Y33MeV_239.dat'
    list = [foil1, foil2]
    lambda_ = np.log(2)/(106.626*d)
    return list, lambda_

def Y_24Na(): #non-mon, single
    foil1 = path + '24Na_Y33MeV_239.dat'
    list = [foil1]
    lambda_ = np.log(2)/(14.997*h)
    return list, lambda_



### INDIUM FOILS ###

def In_116mIn(): #non-mon, single
    foil1 = path + '116mIn_In16MeV_149.dat'
    foil2 = path + '116mIn_In33MeV_349.dat'
    list = [foil1, foil2]
    lambda_ = np.log(2)/(54.29*m)
    return list, lambda_


def In_112mIn(): #non-mon, single
    foil1 = path + '112mIn_In33MeV_349.dat'
    list = [foil1]
    lambda_ = np.log(2)/(20.67*m)
    return list, lambda_


def In_112In(): #double decay # LEGG INN CSV FILEN
    foil1 = path + '112In_In16MeV_349.dat'
    foil2 = path + '112In_In33MeV_149.dat'
    list = [foil1, foil2]
    lambda_parent = np.log(2)/(20.67*m)
    lambda_daughter = np.log(2)/(49.51*d)
    return list, lambda_parent, lambda_daughter


def In_111In(): #non-mon, single
    foil1 = path + '111In_In33MeV_349.dat'
    list = [foil1]
    lambda_ = np.log(2)/(2.8*d)
    return list, lambda_


def In_115mIn(): #non-mon, single
    foil1 = path + '115mIn_In16MeV_149.dat'
    foil2 = path + '115mIn_In33MeV_349.dat'
    list = [foil1, foil2]
    lambda_ = np.log(2)/(4.486*h)
    return list, lambda_


def In_113mIn(): #non-mon, single
    foil1 = path + '113mIn_In16MeV_149.dat'
    foil2 = path + '113mIn_In33MeV_349.dat'
    list = [foil1, foil2]
    lambda_ = np.log(2)/(99.476*m)
    return list, lambda_


def In_114mIn(): #non-mon, single
    foil1 = path + '114mIn_In16MeV_149.dat'
    foil2 = path + '114mIn_In33MeV_349.dat'
    list = [foil1, foil2]
    lambda_ = np.log(2)/(49.51*d)
    return list, lambda_


def In_24Na(): #non-mon, single
    foil1 = path + '24Na_In33MeV_349.dat'
    list = [foil1]
    lambda_ = np.log(2)/(14.997*h)
    return list, lambda_



### ALUMINUM FOILS ###

def Al_24Na(): #non-mon, single
    foil1 = path + '24Na_Al16MeV_113.dat'
    foil2 = path + '24Na_Al33MeV_213.dat'
    list = [foil1, foil2]
    lambda_ = np.log(2)/(14.997*h)
    return list, lambda_


#def Al_(): #non-mon, single
#    foil1 = path + ''
#    foil2 = path + ''
#    list = [foil1, foil2]
#    lambda_ = np.log(2)/(72.1*m)
#    return list, lambda_
