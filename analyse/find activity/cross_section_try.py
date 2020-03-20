
import os
import numpy as np
import matplotlib.pyplot as plt


from foil_info import *

path_ = '/Users/Nora/Documents/UiO/Masteroppgaven/analyse/'


def ALICE(foil, A, Z):
    filename = path + '/../Alice/plot_{}'.format(foil)

    with open(filename) as f:

        begin =   'Ebeam  Z   A   Total    grnd st  isomer1 isomer2  isomer1  isomer2  didl nolevs'
        end   =   'Ebeam =  '
        content_full = f.readlines()
        ind_begin = [line for line in range(len(content_full)) if begin in content_full[line]][0]+2  # only one element but want it as integer
        ind_end   = [line for line in range(len(content_full)) if end   in content_full[line]][0]-1  # list of different, only want the first element
        content = content_full[ind_begin:ind_end]
        #print(content_full[ind_begin])
        #print("ind_begin: ", ind_begin)
        #print("ind_end: ", ind_end)


        E  = np.genfromtxt(filename, delimiter=' ', usecols=[0],skip_header=59, skip_footer=(len(content_full)-len(content)))
        #Z_  = np.genfromtxt(filename, delimiter=' ', usecols=[1], skip_header=56, skip_footer=(len(content_full)-len(content)))
        #A_  = np.genfromtxt(filename, delimiter=' ', usecols=[2], skip_header=56, skip_footer=(len(content_full)-len(content)))
        #print(E)
        CS = np.genfromtxt(filename, delimiter=' ', usecols=[5], skip_header=59, skip_footer=(len(content_full)-len(content)))
        print(CS)
        E_new = []; CS_new = []
        Z = ' ' + Z + ' '
        A = ' ' + A + ' '

        for lines in range(len(content)):
            #print(lines)
            if Z in content[lines] and A in content[lines]:
                #print(Z,A)
                E_new.append(E[lines])
                CS_new.append(CS[lines])
                #print(Z, A)
                #print(content[lines])
                #print("E: ",E_new)
                #print("CS: ", CS_new)
            #else:
            #    print(content[lines])

            #    print(A,Z)
            #    print('---------------------')
                #print('nope')

        #f.close()
        #print("E: ",E_new)
        #print("CS: ", CS_new)
        plt.plot(E_new, CS_new, label='ALICE')
        #plt.show()
        return E, CS


path_ = '/Users/Nora/Documents/UiO/Masteroppgaven/analyse/'


def TALYS(foil, Z, A, file_ending='.tot'):
    # Z = 0XX, A=0XX

    filename = path + '../Talys/' +foil+ '/rp'+Z+A+ file_ending #'.tot'
    #filename = self.path + '/../Talys/' +foil+ '/rp'+Z+A+'.L02'
    E  = np.genfromtxt(filename, delimiter=' ', usecols=[0],skip_header=5)
    CS = np.genfromtxt(filename, delimiter=' ', usecols=[1],skip_header=5)

    #print(CS)
    #print('TALYS CS: ', CS)
    #print(E)

    plt.plot(E,CS, 'r', label='TALYS')

    #plt.show()
    return E, CS


def EXFOR(foil, reaction):
    filename = path + '/../EXFOR/' + foil +'/' + reaction + '.txt'
    #print(filename)
    print(reaction)


    with open(filename) as f:
        begin =   'EXFOR-ID'
        end   =   '//'
        content_full = f.readlines()
        #print(content_full)
        ind_begin = [line for line in range(len(content_full)) if begin in content_full[line]][0]+1  # only one element but want it as integer
        ind_end   = [line for line in range(len(content_full)) if end   in content_full[line]][0]  # list of different, only want the first element
        #print(ind_begin)
        #print(content_full[ind_begin])
        #print(content_full[ind_end])
        #print(ind_end)
        #print(content_full[ind_end])
        content = content_full[ind_begin:ind_end]
        #print(content)

        #str1 = content[0]

        #x = str1.lstrip()
        #print(x.split())
        E = []; dE=[]; CS = []; dCS=[]; author=[]

        for ind in range(len(content)):
            string= content[ind]
            string = (string.lstrip()).split()
            E.append(float(string[0]))
            dE.append(float(string[1]))
            CS.append(float(string[2])*1e3) # in mb
            dCS.append(float(string[3])*1e3) # in mb
            author.append(string[5]) #index 4 is equal to #

        #print('EXFOR CS: ', CS)
        #plt.plot(E, CS)
            #plt.show()
                #plt.errorbar(E[i], CS[i], marker='.', linewidth=0.001, xerr=dE[ind], yerr=dCS[ind], elinewidth=0.5, capthick=0.5, capsize=3.0, label=author[ind] )
        plt.errorbar(E, CS, marker='.', linewidth=0.001, xerr=dE, yerr=dCS, elinewidth=0.5, capthick=0.5, capsize=3.0, label=author[0] )
        plt.legend()
        #plt.show()

        return E, dE, CS, dCS, author
    # else:
    #     print("exfor file does not exist for {}".format(reaction))
    #     return 0, 0, 0, 0, '0'


def Tendl(self, foil, A, Z, file_ending='.tot'):

    #print("foil: ",foil )
    #print("Z: ", Z )
    #print("A: ", A  )

    if foil == 'Ir':
        #A = ['191', '193'] # stable iridium isotopes
        abund_191Ir = 0.373 ; abund_193Ir = 0.627
        #file_ending =
        #endre f_191I til de stabile isotopene i zink osv
        f_191Ir = self.path + '/../Tendl/' + foil + '/rp077191_' + Z+ A + file_ending + '.txt'
        f_193Ir = self.path + '/../Tendl/' + foil + '/rp077193_' + Z +A + file_ending + '.txt'

        #print("Ir 193 file: ",f_193Ir)
        #print("Ir 191 file: ",f_191Ir)
        if os.path.isfile(f_191Ir):
            #print("Ir 191 file: ",f_191Ir)
            #print("f_191Ir exists")
            CS_191Ir = np.genfromtxt(f_191Ir, delimiter=' ', usecols=[1],skip_header=5)
            E_191Ir = np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        else:
            #print("Ir 191 file does not exist")
            CS_191Ir = 0
            E_191Ir =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        if os.path.isfile(f_193Ir):
            #print("f_193Ir exists")
            CS_193Ir = np.genfromtxt(f_193Ir, delimiter=' ', usecols=[1],skip_header=5)
            E_193Ir = np.genfromtxt(f_193Ir, delimiter=' ', usecols=[0],skip_header=5)

        else:
            #print("Ir 193 file does not exist")
            CS_193Ir = 0
            E_193Ir = 0#np.genfromtxt(f_193Ir, delimiter=' ', usecols=[0],skip_header=5)

        E = E_191Ir*abund_191Ir + E_193Ir*abund_193Ir
        CS = CS_191Ir*abund_191Ir + CS_193Ir*abund_193Ir


        #plt.plot(E, CS, label='tendl')




        #plt.plot(E_191Ir,CS_191Ir, label='191Ir')
        #plt.plot(E_193Ir,CS_193Ir, label='193Ir')
        #plt.plot(E, CS, label='tot')
        #plt.legend()
        #plt.show()
        return E, CS

        #except:
            #print("files not exist or file ending is wrong. ")


        #plt.legend()
        #plt.show()

        #E = E_191Ir*abund_191Ir + E_193Ir*abund_193Ir
        #CS = CS_191Ir*abund_191Ir + CS_193Ir*abund_193Ir

    elif foil == 'Ni':
        #finn abundance til hvert isotop fra nndc
        abund_58Ni = 0.68077; abund_60Ni = 0.26233; abund_61Ni = 0.011399; abund_62Ni = 0.036346; abund_64Ni = 0.009255;


        f_58Ni = self.path + '/../Tendl/' + foil + '/rp028058_' + Z + A + file_ending + '.txt'
        f_60Ni = self.path + '/../Tendl/' + foil + '/rp028060_' + Z + A + file_ending + '.txt'
        f_61Ni = self.path + '/../Tendl/' + foil + '/rp028061_' + Z + A + file_ending + '.txt'
        f_62Ni = self.path + '/../Tendl/' + foil + '/rp028062_' + Z + A + file_ending + '.txt'
        f_64Ni = self.path + '/../Tendl/' + foil + '/rp028064_' + Z + A + file_ending + '.txt'

        #print(f_60Ni)
        #print(f_61Ni)
        #print(f_62Ni)
        #print(f_64Ni)
        #E = []; CS = []


        if os.path.isfile(f_58Ni):
            #print("Ir 191 file: ",f_191Ir)
            print("f_58Ni exists")
            CS_58Ni = np.genfromtxt(f_58Ni, delimiter=' ', usecols=[1],skip_header=5)
            E_58Ni = np.genfromtxt(f_58Ni, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("Ni 58 file does not exist")
            CS_58Ni = 0
            E_58Ni =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        if os.path.isfile(f_60Ni):
            #print("Ir 191 file: ",f_191Ir)
            print("f_60Ni exists")
            CS_60Ni = np.genfromtxt(f_60Ni, delimiter=' ', usecols=[1],skip_header=5)
            E_60Ni = np.genfromtxt(f_60Ni, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("Ni 60 file does not exist")
            CS_60Ni = 0
            E_60Ni =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        if os.path.isfile(f_61Ni):
            #print("Ir 191 file: ",f_191Ir)
            print("f_61Ni exists")
            CS_61Ni = np.genfromtxt(f_61Ni, delimiter=' ', usecols=[1],skip_header=5)
            E_61Ni = np.genfromtxt(f_61Ni, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("Ni 61 file does not exist")
            CS_61Ni = 0
            E_61Ni =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        if os.path.isfile(f_62Ni):
            #print("Ir 191 file: ",f_191Ir)
            print("f_62Ni exists")
            CS_62Ni = np.genfromtxt(f_62Ni, delimiter=' ', usecols=[1],skip_header=5)
            E_62Ni = np.genfromtxt(f_62Ni, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("Ni 62 file does not exist")
            CS_62Ni = 0
            E_62Ni =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        if os.path.isfile(f_64Ni):
            #print("Ir 191 file: ",f_191Ir)
            print("f_64Ni exists")
            CS_64Ni = np.genfromtxt(f_64Ni, delimiter=' ', usecols=[1],skip_header=5)
            E_64Ni = np.genfromtxt(f_64Ni, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("Ni 64 file does not exist")
            CS_64Ni = 0
            E_64Ni =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        #CS_58Ni = np.genfromtxt(f_58Ni, delimiter=' ', usecols=[1],skip_header=5)
        #E_58Ni = np.genfromtxt(f_58Ni, delimiter=' ', usecols=[0],skip_header=5)
        #CS_60Ni = np.genfromtxt(f_60Ni, delimiter=' ', usecols=[1],skip_header=5)
        #E_60Ni = np.genfromtxt(f_60Ni, delimiter=' ', usecols=[0],skip_header=5)

        #CS_61Ni = np.genfromtxt(f_61Ni, delimiter=' ', usecols=[1],skip_header=5)
        #E_61Ni = np.genfromtxt(f_61Ni, delimiter=' ', usecols=[0],skip_header=5)
        #CS_62Ni = np.genfromtxt(f_62Ni, delimiter=' ', usecols=[1],skip_header=5)
        #E_62Ni = np.genfromtxt(f_62Ni, delimiter=' ', usecols=[0],skip_header=5)
        #CS_64Ni = np.genfromtxt(f_64Ni, delimiter=' ', usecols=[1],skip_header=5)
        #E_64Ni = np.genfromtxt(f_64Ni, delimiter=' ', usecols=[0],skip_header=5)
        #print(E_64Ni)

        #print(E_58Ni)
        #CS_193Ir = np.genfromtxt(f_193Ir, delimiter=' ', usecols=[1],skip_header=5)
        #E_193Ir = np.genfromtxt(f_193Ir, delimiter=' ', usecols=[0],skip_header=5)

        E = E_58Ni*abund_58Ni + E_60Ni*abund_60Ni + E_61Ni*abund_61Ni + E_62Ni*abund_62Ni + E_64Ni*abund_64Ni
        CS = CS_58Ni*abund_58Ni + CS_60Ni*abund_60Ni + CS_61Ni*abund_61Ni + CS_62Ni*abund_62Ni + CS_64Ni*abund_64Ni
        #CS = CS_191Ir*abund_191Ir + CS_193Ir*abund_193Ir

        #plt.plot(E_60Ni,CS_60Ni, label='60Ni', linewidth=0.5)
        #plt.plot(E_61Ni,CS_61Ni, label='61Ni', linewidth=0.5)
        #plt.plot(E_62Ni,CS_62Ni, label='62Ni', linewidth=0.5)
        #plt.plot(E_64Ni,CS_64Ni, label='64Ni', linewidth=0.5)
        #plt.plot(E, CS, label='tot')

        #except:

            #print("files not exist or file ending is wrong. ")
            #pass

        #plt.legend()
        #plt.show()

    return E, CS






def mydata(filename_CS, foil):
    import os
    print(os.getcwd())
    print('****************')
    file =  './Cross_sections_csv/' + foil + '/' + filename_CS  #'.tot'
    print(file)

    with open(file) as f:
        begin = f.readlines()[1:]
        CS_energi = []; CS_zn_64Cu = []
        for line in begin:
            lines = line.split(',')
            print(lines)
            print('xxxxxxxxxxxxxxxxxxx')
            CS_energi.append(float(lines[0]))
            CS_zn_64Cu.append(float(lines[1]))
            #print('E :', CS_energi)
            #print('CS :', CS_zn_64Cu)

            #for i in range(len(begin)):
            #    string= begin[i]
            #    string = (string.lstrip()).split()
            #    print(string)
            #    print('<<<<<<<<<<<<<<<<')
            #    CS_energi.append(string[0])
            #    CS_zn_64Cu.append(string[1])

    #print(CS_energi)
        print('E :', CS_energi)
        print('CS :', CS_zn_64Cu)
        plt.plot(CS_energi, CS_zn_64Cu,'m.', label='Mydata')

    #plt.show()

    return CS_energi, CS_zn_64Cu


def Cross_section(foil, A, Z, reaction, filename_CS, file_ending='.tot'):

    ALICE(foil, A, Z)
    if len(A) == 2:
        A = '0' + A
    if len(Z) == 2:
        Z = '0' + Z
    TALYS(foil, Z, A)
    EXFOR(foil, reaction)
    mydata(filename_CS, foil)

    plt.legend()
    plt.title(reaction)
    plt.savefig('Cross_section_curves/'+ reaction + '.png', dpi=300)
    plt.show()

### ZINK ###

#Cross_section('Zn', '62', '30', 'Zn(n,x)62Zn', '62ZN')
#Cross_section('Zn', '63', '30', 'Zn(n,x)63Zn', '63ZN')
#Cross_section('Zn', '64', '29', 'Zn(n,x)64Cu', '64CU')
#Cross_section('Zn', '65', '28', 'Zn(n,x)65Ni', '65NI')
#Cross_section('Zn', '65', '30', 'Zn(n,x)65Zn', '65ZN')
#Cross_section('Zn', '66', '29', 'Zn(n,x)66Cu', '66CU')
#Cross_section('Zn', '66', '28', '70Zn(n,na)66Ni', '66NI')
#Cross_section('Zn', '67', '29', '67Zn(n,p)Cu67', '67CU')
Cross_section('Zn', '69', '30', '70Zn(n,2n)69mZn', '69ZNm')






# ALICE('Zn', '64', '29')
# TALYS('Zn', '029', '064')
# EXFOR('Zn(n,x)64Cu')
# mydata('zn_64Cu')

#plt.legend()
#plt.show()
