label='TENDL'
import os
import numpy as np
import matplotlib.pyplot as plt
#import sympy  as sy
#import integrate




from foil_info import *

path_ = '/Users/Nora/Documents/UiO/Masteroppgaven/analyse/'

#xaxis = np.linspace(0,40,1000)

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
        #print(CS)
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
        #plt.plot(E_new, xaxis, label='ALICE')
        plt.plot(E_new, CS_new, color='y', label='ALICE')
        #plt.show()
        return E, CS


path_ = '/Users/Nora/Documents/UiO/Masteroppgaven/analyse/'


def TALYS(foil, Z, A, file_ending='.tot'):
    # Z = 0XX, A=0XX

    #Zn
    filename = path + '../Talys/' +foil+ '/rp'+ Z + A + file_ending #'.tot'
    #zr
    #filename = path + '../Talys/' +foil+ '/rp0'+ Z + '0' + A + file_ending #'.tot'

    #filename = self.path + '/../Talys/' +foil+ '/rp'+Z+A+'.L02'
    E  = np.genfromtxt(filename, delimiter=' ', usecols=[0],skip_header=5)
    CS = np.genfromtxt(filename, delimiter=' ', usecols=[1],skip_header=5)

    #print(CS)
    #print('TALYS CS: ', CS)
    #print(E)

    plt.plot(E, CS, color='g', label='TALYS')

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
            #plt.errorbar(E[ind], CS[ind], marker='.', linewidth=0.001, xerr=dE[ind], yerr=dCS[ind], elinewidth=0.5, capthick=0.5, capsize=3.0, label=author[ind] )
        plt.errorbar(E, CS, marker='.', linewidth=0.001, xerr=dE, yerr=dCS, elinewidth=0.5, capthick=0.5, capsize=3.0, label=author[0] )
        #plt.legend()
        #plt.show()

        return E, dE, CS, dCS, author
    # else:
    #     print("exfor file does not exist for {}".format(reaction))
    #     return 0, 0, 0, 0, '0'


def Tendl(foil, A, Z, file_ending='.tot'):

    #print("foil: ",foil )
    #print("Z: ", Z )
    #print("A: ", A  )

    if foil == 'Zn':
        #A = ['191', '193'] # stable iridium isotopes
        abund_64Zn = 0.4917 ; abund_66Zn = 0.2773 ; abund_67Zn = 0.404 ; abund_68Zn = 0.1845 ; abund_70Zn = 0.061;
        #file_ending =
        #endre f_191I til de stabile isotopene i zink osv
        f_64Zn = path + '/../Tendl/' + foil + '/rp030064_' + Z+ A + file_ending + '.txt'
        f_66Zn = path + '/../Tendl/' + foil + '/rp030066_' + Z +A + file_ending + '.txt'
        f_67Zn = path + '/../Tendl/' + foil + '/rp030067_' + Z +A + file_ending + '.txt'
        f_68Zn = path + '/../Tendl/' + foil + '/rp030068_' + Z +A + file_ending + '.txt'
        f_70Zn = path + '/../Tendl/' + foil + '/rp030070_' + Z +A + file_ending + '.txt'


        #print("Ir 193 file: ",f_193Ir)
        #print("Ir 191 file: ",f_191Ir)
        if os.path.isfile(f_64Zn):
            #print("Ir 191 file: ",f_191Ir)
            #print("f_191Ir exists")
            CS_64Zn = np.genfromtxt(f_64Zn, delimiter=' ', usecols=[1],skip_header=5)
            E_64Zn = np.genfromtxt(f_64Zn, delimiter=' ', usecols=[0],skip_header=5)
        else:
            #print("Ir 191 file does not exist")
            CS_64Zn = 0
            E_64Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        if os.path.isfile(f_66Zn):
            #print("Ir 191 file: ",f_191Ir)
            #print("f_191Ir exists")
            CS_66Zn = np.genfromtxt(f_66Zn, delimiter=' ', usecols=[1],skip_header=5)
            E_66Zn = np.genfromtxt(f_66Zn, delimiter=' ', usecols=[0],skip_header=5)
        else:
            #print("Ir 191 file does not exist")
            CS_66Zn = 0
            E_66Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        if os.path.isfile(f_67Zn):
            #print("Ir 191 file: ",f_191Ir)
            #print("f_191Ir exists")
            CS_67Zn = np.genfromtxt(f_67Zn, delimiter=' ', usecols=[1],skip_header=5)
            E_67Zn = np.genfromtxt(f_67Zn, delimiter=' ', usecols=[0],skip_header=5)
        else:
            #print("Ir 191 file does not exist")
            CS_67Zn = 0
            E_67Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        if os.path.isfile(f_68Zn):
            #print("Ir 191 file: ",f_191Ir)
            #print("f_191Ir exists")
            CS_68Zn = np.genfromtxt(f_68Zn, delimiter=' ', usecols=[1],skip_header=5)
            E_68Zn = np.genfromtxt(f_68Zn, delimiter=' ', usecols=[0],skip_header=5)
        else:
            #print("Ir 191 file does not exist")
            CS_68Zn = 0
            E_68Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        if os.path.isfile(f_70Zn):
            #print("Ir 191 file: ",f_191Ir)
            #print("f_191Ir exists")
            CS_70Zn = np.genfromtxt(f_70Zn, delimiter=' ', usecols=[1],skip_header=5)
            E_70Zn = np.genfromtxt(f_70Zn, delimiter=' ', usecols=[0],skip_header=5)
        else:
            #print("Ir 191 file does not exist")
            CS_70Zn = 0
            E_70Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        E = E_64Zn*abund_64Zn + E_66Zn*abund_66Zn + E_67Zn*abund_67Zn + E_68Zn*abund_68Zn + E_70Zn*abund_70Zn
        CS = CS_64Zn*abund_64Zn + CS_66Zn*abund_66Zn + CS_67Zn*abund_67Zn + CS_68Zn*abund_68Zn + CS_70Zn*abund_70Zn

        print(len(E))
        print(len(CS))
        print('...........................')
        #plt.plot(E, xaxis, 'y--', label='TENDL')
        plt.plot(E, CS, color='r', label='TENDL')

        #plt.plot(E_191Ir,CS_191Ir, label='191Ir')
        #plt.plot(E_193Ir,CS_193Ir, label='193Ir')
        #plt.plot(E, CS, label='tot')
        #plt.legend()
        #plt.show()
        return E, CS



    elif foil == 'Zr':
        #finn abundance til hvert isotop fra nndc
        abund_90Zr = 0.5145 ; abund_91Zr = 0.1122 ; abund_92Zr = 0.1715 ; abund_94Zr = 0.1738 ; abund_96Zr = 0.280 ;


        #f_90Zr = path + '/../Tendl/' + foil + '/rp040090_' + Z  + A + file_ending + '.txt'
        f_90Zr = path + '/../Tendl/' + foil + '/rp040090_0' + Z + '0' + A + '.L02.ave' + '.txt'
        f_91Zr = path + '/../Tendl/' + foil + '/rp040091_' + Z  + A + file_ending + '.txt'
        f_92Zr = path + '/../Tendl/' + foil + '/rp040092_' + Z  + A + file_ending + '.txt'
        f_94Zr = path + '/../Tendl/' + foil + '/rp040094_' + Z  + A + file_ending + '.txt'
        f_96Zr = path + '/../Tendl/' + foil + '/rp040096_' + Z  + A + file_ending + '.txt'


        print(f_90Zr)
        #print(f_61Ni)
        #print(f_62Ni)
        #print(f_64Ni)
        #E = []; CS = []


        if os.path.isfile(f_90Zr):
            #print("Ir 191 file: ",f_191Ir)
            print("f_90Zr exists")
            CS_90Zr = np.genfromtxt(f_90Zr, delimiter=' ', usecols=[1],skip_header=5) * abund_90Zr
            print(len(CS_90Zr))
            #E_90Zr = np.genfromtxt(f_90Zr, delimiter=' ', usecols=[0],skip_header=5)
            E_90Zr = np.genfromtxt(f_90Zr, delimiter=' ', usecols=[0],skip_header=5)
            print(len(E_90Zr))
        else:
            print("90Zr file does not exist")
            CS_90Zr = 0
            E_90Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        if os.path.isfile(f_91Zr):
            #print("Ir 191 file: ",f_191Ir)
            print("f_91Zr exists")
            CS_91Zr = np.genfromtxt(f_91Zr, delimiter=' ', usecols=[1],skip_header=5) * abund_91Zr
            print(len(CS_91Zr))
            E_91Zr = np.genfromtxt(f_91Zr, delimiter=' ', usecols=[0],skip_header=5)
            print(len(E_91Zr))
        else:
            print("91Zr file does not exist")
            CS_91Zr = 0
            E_91Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        if os.path.isfile(f_92Zr):
            #print("Ir 191 file: ",f_191Ir)
            print("f_92Zr exists")
            CS_92Zr = np.genfromtxt(f_92Zr, delimiter=' ', usecols=[1],skip_header=5) * abund_92Zr
            print(len(CS_92Zr))
            E_92Zr = np.genfromtxt(f_92Zr, delimiter=' ', usecols=[0],skip_header=5)
            print(len(E_92Zr))
        else:
            print("92Zr file does not exist")
            CS_92Zr = 0
            E_92Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        if os.path.isfile(f_94Zr):
            #print("Ir 191 file: ",f_191Ir)
            print("f_94Zr exists")
            CS_94Zr = np.genfromtxt(f_94Zr, delimiter=' ', usecols=[1],skip_header=5) * abund_94Zr
            print(len(CS_94Zr))
            E_94Zr = np.genfromtxt(f_94Zr, delimiter=' ', usecols=[0],skip_header=5)
            print(len(E_94Zr))

        else:
            print("94Zr file does not exist")
            CS_94Zr = 0
            E_94Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        if os.path.isfile(f_96Zr):
            #print("Ir 191 file: ",f_191Ir)
            print("f_96Zr exists")
            CS_96Zr = np.genfromtxt(f_96Zr, delimiter=' ', usecols=[1],skip_header=5) * abund_96Zr
            print(len(CS_96Zr))
            E_96Zr = np.genfromtxt(f_96Zr, delimiter=' ', usecols=[0],skip_header=5)
            print(len(E_96Zr))

        else:
            print("96Zr file does not exist")
            CS_96Zr = 0
            E_96Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        #print(E_58Ni)
        #CS_193Ir = np.genfromtxt(f_193Ir, delimiter=' ', usecols=[1],skip_header=5)
        #E_193Ir = np.genfromtxt(f_193Ir, delimiter=' ', usecols=[0],skip_header=5)

        E = E_90Zr*abund_90Zr + E_91Zr*abund_91Zr + E_92Zr*abund_92Zr + E_94Zr*abund_94Zr + E_96Zr*abund_96Zr
        #CS = CS_90Zr*abund_90Zr + CS_91Zr*abund_91Zr + CS_92Zr*abund_92Zr + CS_94Zr*abund_94Zr + CS_96Zr*abund_96Zr
        CS = CS_90Zr + CS_91Zr + CS_92Zr + CS_94Zr + CS_96Zr

        print('E_zr :', len(E))
        print('CS_zr :', len(CS))
        #CS = CS_191Ir*abund_191Ir + CS_193Ir*abund_193Ir

        #plt.plot(E, xaxis, 'y--', label='TENDL')

        plt.plot(E, CS, color='r', label='TENDL')
        #plt.plot(E_61Ni,CS_61Ni, label='61Ni', linewidth=0.5)

        return E, CS



    elif foil == 'Y':
        #finn abundance til hvert isotop fra nndc
        abund_89Y = 1.00;


        f_89Y = path + '/../Tendl/' + foil + '/rp_039089' + Z + A + file_ending + '.txt'


        if os.path.isfile(f_89Y):
            #print("Ir 191 file: ",f_191Ir)
            print("f_89Y exists")
            CS_89Y = np.genfromtxt(f_89Y, delimiter=' ', usecols=[1],skip_header=5)
            E_89Y = np.genfromtxt(f_89Y, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("89Y file does not exist")
            CS_90Zr = 0
            E_90Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        E = E_89Y*abund_89Y
        CS = CS_89Y*abund_89Y

        plt.plot(E, CS, color='r', label='TENDL')

        return E, CS


    elif foil == 'In':
        #finn abundance til hvert isotop fra nndc
        abund_113In = 0.429;  abund_115In = 0.9571;


        f_abund_113In = path + '/../Tendl/' + foil + '/rp_049113' + Z + A + file_ending + '.txt'
        f_abund_115In = path + '/../Tendl/' + foil + '/rp_049115' + Z + A + file_ending + '.txt'


        if os.path.isfile(f_113Y):
            #print("Ir 191 file: ",f_191Ir)
            print("f_113Y exists")
            CS_113Y = np.genfromtxt(f_113Y, delimiter=' ', usecols=[1],skip_header=5)
            E_113Y = np.genfromtxt(f_113Y, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("113Y file does not exist")
            CS_113Y = 0
            E_113Y =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        if os.path.isfile(f_115Y):
            #print("Ir 191 file: ",f_191Ir)
            print("f_115Y exists")
            CS_115Y = np.genfromtxt(f_115Y, delimiter=' ', usecols=[1],skip_header=5)
            E_115Y = np.genfromtxt(f_115Y, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("115Y file does not exist")
            CS_115Y = 0
            E_115Y =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        E = E_113Y*abund_113Y + E_115Y*abund_115Y
        CS = CS_113Y*abund_113Y + CS_115Y*abund_115Y

        plt.plot(E, CS, color='r', label='TENDL')

        return E, CS


    elif foil == 'Al':
        #finn abundance til hvert isotop fra nndc
        abund_27Al = 1.00;


        f_abund_27Al = path + '/../Tendl/' + foil + '/rp_013027' + Z + A + file_ending + '.txt'


        if os.path.isfile(f_27Al):
            #print("Ir 191 file: ",f_191Ir)
            print("f_27Al exists")
            CS_27Al = np.genfromtxt(f_27Al, delimiter=' ', usecols=[1],skip_header=5)
            E_27Al = np.genfromtxt(f_27Al, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("27Al file does not exist")
            CS_27Al = 0
            E_27Al =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        E = E_27Al*abund_27Al
        CS = CS_27Al*abund_27Al

        #plt.plot(E, xaxis, 'y--', label='TENDL')

        plt.plot(E, CS, color='r', label='TENDL')

        return E, CS





def CoH(foil, A, Z, filename_CS):

    #print("foil: ",foil )
    #print("Z: ", Z )
    #print("A: ", A  )

    if foil == 'Zn':
        #A = ['191', '193'] # stable iridium isotopes
        abund_64Zn = 0.4917 ; abund_66Zn = 0.2773 ; abund_67Zn = 0.404 ; abund_68Zn = 0.1845 ; abund_70Zn = 0.061;
        #file_ending =
        #endre f_191I til de stabile isotopene i zink osv
        v_64Zn = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-0' + filename_CS + '_coh' + '.txt'
        #v_64Zn_emire = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-' + 'Cu-64' + '_empire' + '.txt'
        v_67Zn = path + '/../EMPIRECOH2/' + foil + '/67Zn/' +  Z + '-0' + filename_CS + '_coh' + '.txt'        #f_67Zn = path + '/../EMPIRECOH2/' + foil + '/rp030067_' + Z +A + file_ending + '.txt'
        #f_68Zn = path + '/../EMPIRECOH2/' + foil + '/rp030068_' + Z +A + file_ending + '.txt'
        #f_70Zn = path + '/../EMPIRECOH2/' + foil + '/rp030070_' + Z +A + file_ending + '.txt'


        #print("Zn 64 file: ",v_64Zn)
        #print("Ir 191 file: ",f_191Ir)
        if os.path.isfile(v_64Zn):
            #print("Ir 191 file: ",f_191Ir)
            print("v_64Zn exists")
            CS_64Zn_Coh = np.genfromtxt(v_64Zn, delimiter='	', usecols=[1])
            E_64Zn_Coh = np.genfromtxt(v_64Zn, delimiter='	', usecols=[0])
        else:
            print("64 Zn file does not exist")
            CS_64Zn_Coh = 0
            E_64Zn_Coh =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        #print('!!!!!!!!!!!!!!!!!!!!!!!!!')
        #print(E_64Zn_Coh)
        plt.plot(E_64Zn_Coh, CS_64Zn_Coh, color='c', label='CoH')


        if os.path.isfile(v_67Zn):
            #print("Ir 191 file: ",f_191Ir)
            print("v_64Zn exists")
            CS_67Zn_Coh = np.genfromtxt(v_67Zn, delimiter='	', usecols=[1])
            E_67Zn_Coh = np.genfromtxt(v_67Zn, delimiter='	', usecols=[0])
        else:
            print("64 Zn file does not exist")
            CS_67Zn_Coh = 0
            E_67Zn_Coh =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        #print('!!!!!!!!!!!!!!!!!!!!!!!!!')
        #print(E_67Zn_Coh)
        #plt.plot(E_67Zn_Coh, CS_67Zn_Coh, color='c', label='CoH')


def EMPIRE(foil, A, Z, filename_CS):

    #print("foil: ",foil )
    #print("Z: ", Z )
    #print("A: ", A  )

    if foil == 'Zn':
        #A = ['191', '193'] # stable iridium isotopes
        abund_64Zn = 0.4917 ; abund_66Zn = 0.2773 ; abund_67Zn = 0.404 ; abund_68Zn = 0.1845 ; abund_70Zn = 0.061;
        #file_ending =
        #endre f_191I til de stabile isotopene i zink osv
        #v_64Zn = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-0' + filename_CS + '_coh' + '.txt'
        u_64Zn_emire = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-' + 'Cu-64' + '_empire' + '.txt'
        u_68Zn_emire = path + '/../EMPIRECOH2/' + foil + '/68Zn/' +  Z + '-' + 'Cu-67' + '_empire' + '.txt'

        #f_66Zn = path + '/../EMPIRECOH2/' + foil + '/rp030066_' + Z +A + '_empire' + '.txt'
        #f_67Zn = path + '/../EMPIRECOH2/' + foil + '/rp030067_' + Z +A + file_ending + '.txt'
        #f_68Zn = path + '/../EMPIRECOH2/' + foil + '/rp030068_' + Z +A + file_ending + '.txt'
        #f_70Zn = path + '/../EMPIRECOH2/' + foil + '/rp030070_' + Z +A + file_ending + '.txt'


        print("Zn 64 file: ",u_64Zn_emire)
        #print("Ir 191 file: ",f_191Ir)
        if os.path.isfile(u_64Zn_emire):
            #print("Ir 191 file: ",f_191Ir)
            print("v_64Zn exists")
            CS_64Zn_E = np.genfromtxt(u_64Zn_emire, delimiter='	', usecols=[1])
            E_64Zn_E = np.genfromtxt(u_64Zn_emire, delimiter='	', usecols=[0])
        else:
            print("64 Zn file does not exist")
            CS_64Zn_E = 0
            E_64Zn_E =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        #print('!!!!!!!!!!!!!!!!!!!!!!!!!')
        #print(E_64Zn_E)
        plt.plot(E_64Zn_E, CS_64Zn_E, color='m', label='EMPIRE')


        if os.path.isfile(u_68Zn_emire):
            #print("Ir 191 file: ",f_191Ir)
            print("u_8Zn exists")
            CS_68Zn_E = np.genfromtxt(u_68Zn_emire, delimiter='	', usecols=[1])
            E_68Zn_E = np.genfromtxt(u_68Zn_emire, delimiter='	', usecols=[0])
        else:
            print("68 Zn file does not exist")
            CS_68Zn_E = 0
            E_68Zn_E =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        #print('!!!!!!!!!!!!!!!!!!!!!!!!!')
        #print(E_68Zn_E)
        #plt.plot(E_68Zn_E, CS_68Zn_E, color='m', label='EMPIRE')



        # if os.path.isfile(f_66Zn):
        #     #print("Ir 191 file: ",f_191Ir)
        #     #print("f_191Ir exists")
        #     CS_66Zn = np.genfromtxt(f_66Zn, delimiter=' ', usecols=[1],skip_header=5)
        #     E_66Zn = np.genfromtxt(f_66Zn, delimiter=' ', usecols=[0],skip_header=5)
        # else:
        #     #print("Ir 191 file does not exist")
        #     CS_66Zn = 0
        #     E_66Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        #
        # if os.path.isfile(f_67Zn):
        #     #print("Ir 191 file: ",f_191Ir)
        #     #print("f_191Ir exists")
        #     CS_67Zn = np.genfromtxt(f_67Zn, delimiter=' ', usecols=[1],skip_header=5)
        #     E_67Zn = np.genfromtxt(f_67Zn, delimiter=' ', usecols=[0],skip_header=5)
        # else:
        #     #print("Ir 191 file does not exist")
        #     CS_67Zn = 0
        #     E_67Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        #
        #
        # if os.path.isfile(f_68Zn):
        #     #print("Ir 191 file: ",f_191Ir)
        #     #print("f_191Ir exists")
        #     CS_68Zn = np.genfromtxt(f_68Zn, delimiter=' ', usecols=[1],skip_header=5)
        #     E_68Zn = np.genfromtxt(f_68Zn, delimiter=' ', usecols=[0],skip_header=5)
        # else:
        #     #print("Ir 191 file does not exist")
        #     CS_68Zn = 0
        #     E_68Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        #
        #
        # if os.path.isfile(f_70Zn):
        #     #print("Ir 191 file: ",f_191Ir)
        #     #print("f_191Ir exists")
        #     CS_70Zn = np.genfromtxt(f_70Zn, delimiter=' ', usecols=[1],skip_header=5)
        #     E_70Zn = np.genfromtxt(f_70Zn, delimiter=' ', usecols=[0],skip_header=5)
        # else:
        #     #print("Ir 191 file does not exist")
        #     CS_70Zn = 0
        #     E_70Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        #
        # E = E_64Zn*abund_64Zn + E_66Zn*abund_66Zn + E_67Zn*abund_67Zn + E_68Zn*abund_68Zn + E_70Zn*abund_70Zn
        # CS = CS_64Zn*abund_64Zn + CS_66Zn*abund_66Zn + CS_67Zn*abund_67Zn + CS_68Zn*abund_68Zn + CS_70Zn*abund_70Zn
        #
        # print(len(E))
        # print(len(CS))
        #print('...........................')
        #plt.plot(E, xaxis, 'y--', label='TENDL')
        #plt.plot(E, CS, 'y--', label='TENDL')







def mydata(filename_CS, foil):
    import os
    #print(os.getcwd())
    print('****************')
    file =  './Cross_sections_csv/' + foil + '/' + filename_CS  #'.tot'
    file_33 = path + '../../Jon code/meulders_33MeV.csv'
    file_16 = path + '../../Jon code/meulders_16MeV.csv'
    print(file_33)


    with open(file_33) as p:
        begin2 = p.readlines()[1:]
        #print(begin2)
        E_33 = []; dPhidE_33 = []; x_values_33 = [6.753406008583691, 21.685666141732284]#dE_33 = [5.9795033493567011, 8.9527567837918909]
        for line2 in begin2:
            lines2 = line2.split(',')
            #print(lines2)
            E_33.append(float(lines2[0]))
            dPhidE_33.append(float(lines2[1]))
            #print('_:_:_:_:_:_:_:_:_:_:_:_:_:_')
            #print('EEE: ', E)
            #print(dPhidE)
        # gives average of dPhidE np.trapz(y,x
        #dPhidE_33_avarage = np.trapz(dPhidE_33, E_33)
        # Multiplying two lists
        res_list_33 = [E_33[i] * dPhidE_33[i] for i in range(len(E_33))]
        E_33_average = np.trapz(res_list_33, E_33) / np.trapz(dPhidE_33, E_33)
        dEl_33 = E_33_average-x_values_33[0]
        dEr_33 = x_values_33[1]-E_33_average   #left and right uncertainty in energy
        print('E_33MeV_average:', E_33_average)


    with open(file_16) as p:
        begin3 = p.readlines()[1:]
                #print(begin2)
        E_16 = []; dPhidE_16 = []; x_values_16 = [2.5311, 11.76469534883721]#dE_16 = [3.4329052022917228, 5.8006901465454863]
        for line3 in begin3:
            lines3 = line3.split(',')
            #print(lines3)
            E_16.append(float(lines3[0]))
            dPhidE_16.append(float(lines3[1]))
                # gives average of dPhidE np.trapz(y,x
                #dPhidE_33_avarage = np.trapz(dPhidE_33, E_33)
        res_list_16 = [E_16[i] * dPhidE_16[i] for i in range(len(E_16))]
        E_16_average = np.trapz(res_list_16, E_16) / np.trapz(dPhidE_16, E_16)
        dEl_16 = E_16_average-x_values_16[0]
        dEr_16 = x_values_16[1]-E_16_average
        print('E_16MeV_average:', E_16_average)


        #print(E[i])
        #print(dE_)
        #print('dE:', dE)
        #print(len(dE), dE)
        CS_energi = [E_16_average, E_33_average]
        print(CS_energi)


    with open(file) as f:
        begin = f.readlines()[1:]
        #CS_energi = [];
        CS_zn = []; dCS = []; dE = []
        for line in begin:
            lines = line.split(',')
            #print(lines)
            #print('xxxxxxxxxxxxxxxxxxx')
            #CS_energi.append(float(lines[0]))
            CS_zn.append(float(lines[1]))
            dCS.append(float(lines[2]))
            #dPhidE =
            #dE = int(CS_energi*dPhidE) / int(dPhidE)
            #print('dE :', dE)
            #print('--------------')
            #print('CS :', CS_zn_64Cu)
            dE = np.array(([dEl_16, dEr_16], [dEl_33, dEr_33]))

    #print(CS_energi)
        print('E :', CS_energi)
        print('CS :', CS_zn)
        print('dE :', dE)
        #plt.plot(CS_energi, xaxis, 'g.', label='THIS WORK')
        #plt.plot(CS_energi, CS_zn, 'k.', label='This Work')
        #plt.errorbar(CS_energi, CS_zn, marker='P', color='darkred',linewidth=0.0001, yerr=dCS, elinewidth=1.0, capthick=1.0, capsize=3.0, label='This Work')
        plt.errorbar(CS_energi, CS_zn, marker='.', linewidth=0.001, xerr=dE, yerr=dCS, elinewidth=0.5, capthick=0.5, capsize=3.0, label='This Work' )
        #plt.errorbar(CS_energi[1], CS_zn[1], marker='.', linewidth=0.001, xerr=dE[1], yerr=dCS[1], elinewidth=0.5, capthick=0.5, capsize=3.0, label='This Work' )
        #plt.errorbar(CS_energi, CS_zn, marker='.', linewidth=0.001, xerr=dE, yerr=dCS, elinewidth=0.5, capthick=0.5, capsize=3.0 )

    #plt.show()

    return CS_energi, CS_zn, dCS, dE


def Cross_section(foil, A, Z, reaction, filename_CS, file_ending='.tot'):

    ALICE(foil, A, Z)
    if len(A) == 2:
       A = '0' + A
    if len(Z) == 2:
       Z = '0' + Z
    TALYS(foil, Z, A)

    EXFOR(foil, reaction)
    #no EXFOR for :
    # natZn(n,x)67Cu
    Tendl(foil, A, Z)
    CoH(foil, A, Z, filename_CS)
    EMPIRE(foil, A, Z, filename_CS)

    mydata(filename_CS, foil)

    plt.legend()
    plt.title(reaction)
    plt.xlim(0,40)
    plt.ylim(ymin=0)
    plt.xlabel('Neutron Energy (MeV)')
    plt.ylabel('Cross Section (mb)')
    plt.savefig('Cross_section_curves/'+ reaction + '.png', dpi=300)
    plt.show()


#Cross_section('foil', 'A', 'Z', 'EXFOR', 'mydata')
#Cross_section('foil', 'A', 'Z', 'reaction', 'filename_CS')

### ZINK ###

#Cross_section('Zn', '62', '30', 'Zn(n,x)62Zn', '62ZN')
#Cross_section('Zn', '63', '30', 'Zn(n,x)63Zn', '63ZN')
Cross_section('Zn', '64', '29', 'nat-Zn(n,x)64Cu', '64CU')
#Cross_section('Zn', '65', '28', 'Zn(n,x)65Ni', '65NI')
#Cross_section('Zn', '65', '30', 'Zn(n,x)65Zn', '65ZN')
#Cross_section('Zn', '66', '29', 'Zn(n,x)66Cu', '66CU')
#Cross_section('Zn', '67', '29', 'nat-Zn(n,x)67Cu', '67CU') # really 68Zn not naturall in exfor
#Cross_section('Zn', '66', '28', '70Zn(n,na)66Ni', '66NI')
#Cross_section('Zn', '67', '29', '67Zn(n,p)Cu67', '67CU')
#Cross_section('Zn', '69', '30', '70Zn(n,2n)69mZn', '69ZNm')


### Zirconium ###

#Cross_section('Zr', '90', '39', '90Zr(n,x)90mY', '90Ym')
#Cross_section('Zr', '90', '39', '91Zr(n,np)90mY', '90Ym')
#Cross_section('Zr', '91', '38', '94Zr(n,a)91Sr', '91SR')
#Cross_section('Zr', '91', '39', '91Zr(n,p)91mY', '91Ym')
#Cross_section('Zr', '91', '39', '91zr(n,x)91mY', '91Ym')
#Cross_section('Zr', '92', '39', '92Zr(n,p)92Y', '92Y')
#Cross_section('Zr', '93', '39', '94Zr(n,np)93Y', '93Y')
#Cross_section('Zr', '93', '39', '94Zr(n,x)93Y', '93Y')
#Cross_section('Zr', '95', '40', '94Zr(n,g)95Zr', '95Zr')
#Cross_section('Zr', '95', '40', '96Zr(n,2n)95Zr', '95Zr')
#Cross_section('Zr', '97', '40', '96Zr(n,g)97Zr', '97Zr')

#Cross_section('Zr', '93', '39', '94Zr(n,x)93Y', '95NB') INGEN EXFOR ELLER TENDL
#Cross_section('Zr', '93', '39', '94Zr(n,x)93Y', '97NB') INGEN EXFOR ELLER TENDL
#Cross_section('Zr', '98', '40', '96Zr(n,g)97Zr', '98Zr') INGEN EXFOR ELLER TENDL


### Indium ###

#Cross_section('In', '62', '30', 'Zn(n,x)62Zn', '62ZN')

###  Aluminum ###

#Cross_section('Al', '62', '30', 'Zn(n,x)62Zn', '62ZN')

### Yttrium ###

#Cross_section('Y', '62', '30', 'Zn(n,x)62Zn', '62ZN')



# ALICE('Zn', '64', '29')
# TALYS('Zn', '029', '064')
# EXFOR('Zn(n,x)64Cu')
# mydata('zn_64Cu')

#plt.legend()
#plt.show()
