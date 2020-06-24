label='TENDL'
import os
import numpy as np
import matplotlib.pyplot as plt
#import sympy  as sy
#import integrate

# LINE 22 - ALICE
# LINE 77 - TALYS
# LINE 100 - EXFOR
# LINE 150 - TENDL
# LINE 441 - CoH
# LINE 489 - EMPIRE


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


def TALYS(foil, A, Z, file_ending='.tot'):
    # Z = 0XX, A=0XX

    #Zn
    #filename = path + '../Talys/' +foil+ '/rp'+ Z + A + file_ending #'.tot'
    #Al
    #filename = path + '../Talys/' +foil+ '/rp'+ Z + A + file_ending + '.027' #'.tot'
    # In
    #filename = path + '../Talys/' +foil+ '/rp'+ Z + A + file_ending + '.113' #'.tot
    # Y
    #filename = path + '../Talys/' +foil+ '/rp'+ Z + A + file_ending + '.089' #'.tot
    # Zr
    filename = path + '../Talys/' +foil+ '/rp'+ Z + A + file_ending #'.tot



    print('TALYS filename :', filename)

    #filename = self.path + '/../Talys/' +foil+ '/rp'+Z+A+'.L02'
    E  = np.genfromtxt(filename, delimiter=' ', usecols=[0],skip_header=5)
    CS = np.genfromtxt(filename, delimiter=' ', usecols=[1],skip_header=5)

    #print(CS)
    print('TALYS CS: ', CS)
    print(E)

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
        f_64Zn = path + '/../Tendl/' + foil + '/rp030064_' + Z + A + file_ending + '.txt'
        f_66Zn = path + '/../Tendl/' + foil + '/rp030066_' + Z + A + file_ending + '.txt'
        f_67Zn = path + '/../Tendl/' + foil + '/rp030067_' + Z + A + file_ending + '.txt'
        f_68Zn = path + '/../Tendl/' + foil + '/rp030068_' + Z + A + file_ending + '.txt'
        f_70Zn = path + '/../Tendl/' + foil + '/rp030070_' + Z + A + file_ending + '.txt'


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


        #print(f_90Zr)
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
            CS_89Y = 0
            E_89Y =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        E = E_89Y*abund_89Y
        CS = CS_89Y*abund_89Y

        plt.plot(E, CS, color='r', label='TENDL')

        return E, CS


    elif foil == 'In':
        #finn abundance til hvert isotop fra nndc
        abund_113In = 0.429;  abund_115In = 0.9571;


        f_113In = path + '/../Tendl/' + foil + '/rp_049113_' + Z + A + file_ending + '.txt'
        f_115In = path + '/../Tendl/' + foil + '/rp_049115_' + Z + A + file_ending + '.txt'


        if os.path.isfile(f_113In):
            #print("Ir 191 file: ",f_191Ir)
            print("f_113In exists")
            CS_113In = np.genfromtxt(f_113In, delimiter=' ', usecols=[1],skip_header=5)
            E_113In = np.genfromtxt(f_113In, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("f_113In file does not exist")
            CS_113In = 0
            E_113In =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        if os.path.isfile(f_115In):
            #print("Ir 191 file: ",f_191Ir)
            print("f_115In exists")
            CS_115In = np.genfromtxt(f_115In, delimiter=' ', usecols=[1],skip_header=5)
            E_115In = np.genfromtxt(f_115In, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("f_115In file does not exist")
            CS_115In = 0
            E_115In =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        E = E_113In*abund_113In + E_115In*abund_115In
        CS = CS_113In*abund_113In + CS_115In*abund_115In

        plt.plot(E, CS, color='r', label='TENDL')

        return E, CS


    elif foil == 'Al':
        #finn abundance til hvert isotop fra nndc
        abund_27Al = 1.00;


        f_27Al = path + '/../Tendl/' + foil + '/rp_013027' + Z + A + file_ending + '.txt'


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





# def CoH(foil, A, Z, filename_CS):
#
#     #print("foil: ",foil )
#     #print("Z: ", Z )
#     #print("A: ", A  )
#
#     if foil == 'Zn':
#         #A = ['191', '193'] # stable iridium isotopes
#         abund_64Zn = 0.4917 ; abund_66Zn = 0.2773 ; abund_67Zn = 0.404 ; abund_68Zn = 0.1845 ; abund_70Zn = 0.061;
#         #file_ending =
#         #endre f_191I til de stabile isotopene i zink osv
#         v_64Zn = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-0' + filename_CS + '_coh' + '.txt'
#         #v_64Zn_empire = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-' + 'Cu-64' + '_empire' + '.txt'
#         v_67Zn = path + '/../EMPIRECOH2/' + foil + '/67Zn/' +  Z + '-0' + filename_CS + '_coh' + '.txt'        #f_67Zn = path + '/../EMPIRECOH2/' + foil + '/rp030067_' + Z +A + file_ending + '.txt'
#         #f_68Zn = path + '/../EMPIRECOH2/' + foil + '/rp030068_' + Z +A + file_ending + '.txt'
#         #f_70Zn = path + '/../EMPIRECOH2/' + foil + '/rp030070_' + Z +A + file_ending + '.txt'
#
#
#         #print("Zn 64 file: ",v_64Zn)
#         #print("Ir 191 file: ",f_191Ir)
#         if os.path.isfile(v_64Zn):
#             #print("Ir 191 file: ",f_191Ir)
#             print("v_64Zn exists")
#             CS_64Zn_Coh = np.genfromtxt(v_64Zn, delimiter='	', usecols=[1])
#             E_64Zn_Coh = np.genfromtxt(v_64Zn, delimiter='	', usecols=[0])
#         else:
#             print("64 Zn file does not exist")
#             CS_64Zn_Coh = 0
#             E_64Zn_Coh =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
#         #print('!!!!!!!!!!!!!!!!!!!!!!!!!')
#         #print(E_64Zn_Coh)
#         plt.plot(E_64Zn_Coh, CS_64Zn_Coh, color='c', label='CoH')
#
#
#         if os.path.isfile(v_67Zn):
#             #print("Ir 191 file: ",f_191Ir)
#             print("v_64Zn exists")
#             CS_67Zn_Coh = np.genfromtxt(v_67Zn, delimiter='	', usecols=[1])
#             E_67Zn_Coh = np.genfromtxt(v_67Zn, delimiter='	', usecols=[0])
#         else:
#             print("64 Zn file does not exist")
#             CS_67Zn_Coh = 0
#             E_67Zn_Coh =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
#         #print('!!!!!!!!!!!!!!!!!!!!!!!!!')
#         #print(E_67Zn_Coh)
#         #plt.plot(E_67Zn_Coh, CS_67Zn_Coh, color='c', label='CoH')


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
        # 64Zn
        u_64Zn_empire1 = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-' + 'Cu-64' + '_empire' + '.txt'
        u_64Zn_empire2 = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-' + 'Zn-62' + '_empire' + '.txt'
        u_64Zn_empire3 = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-' + 'Zn-63' + '_empire' + '.txt'
        # 66 Zn
        u_66Zn_empire1 = path + '/../EMPIRECOH2/' + foil + '/66Zn/' +  Z + '-' + 'Cu-64' + '_empire' + '.txt'
        u_66Zn_empire2 = path + '/../EMPIRECOH2/' + foil + '/66Zn/' +  Z + '-' + 'Cu-66' + '_empire' + '.txt'
        u_66Zn_empire3 = path + '/../EMPIRECOH2/' + foil + '/66Zn/' +  Z + '-' + 'Zn-65' + '_empire' + '.txt'
        # 67Zn
        u_67Zn_empire1 = path + '/../EMPIRECOH2/' + foil + '/67Zn/' +  Z + '-' + 'Cu-67' + '_empire' + '.txt'
        u_67Zn_empire2 = path + '/../EMPIRECOH2/' + foil + '/67Zn/' +  Z + '-' + 'Ni-65' + '_empire' + '.txt'
        u_67Zn_empire3 = path + '/../EMPIRECOH2/' + foil + '/67Zn/' +  Z + '-' + 'Ni-66' + '_empire' + '.txt'
        u_67Zn_empire4 = path + '/../EMPIRECOH2/' + foil + '/67Zn/' +  Z + '-' + 'Cu-64' + '_empire' + '.txt'
        u_67Zn_empire5 = path + '/../EMPIRECOH2/' + foil + '/67Zn/' +  Z + '-' + 'Cu-66' + '_empire' + '.txt'
        u_67Zn_empire6 = path + '/../EMPIRECOH2/' + foil + '/67Zn/' +  Z + '-' + 'Cu-67' + '_empire' + '.txt'
        u_67Zn_empire7 = path + '/../EMPIRECOH2/' + foil + '/67Zn/' +  Z + '-' + 'Zn-65' + '_empire' + '.txt'
        # 86Zn
        u_68Zn_empire1 = path + '/../EMPIRECOH2/' + foil + '/68Zn/' +  Z + '-' + 'Ni-65' + '_empire' + '.txt'
        u_68Zn_empire2 = path + '/../EMPIRECOH2/' + foil + '/68Zn/' +  Z + '-' + 'Ni-66' + '_empire' + '.txt'
        u_68Zn_empire3 = path + '/../EMPIRECOH2/' + foil + '/68Zn/' +  Z + '-' + 'Cu-64' + '_empire' + '.txt'
        u_68Zn_empire4 = path + '/../EMPIRECOH2/' + foil + '/68Zn/' +  Z + '-' + 'Cu-66' + '_empire' + '.txt'
        u_68Zn_empire5 = path + '/../EMPIRECOH2/' + foil + '/68Zn/' +  Z + '-' + 'Cu-67' + '_empire' + '.txt'
        u_68Zn_empire6 = path + '/../EMPIRECOH2/' + foil + '/68Zn/' +  Z + '-' + 'Zn-65' + '_empire' + '.txt'
        # 70Zn
        u_70Zn_empire1 = path + '/../EMPIRECOH2/' + foil + '/70Zn/' +  Z + '-' + 'Ni-65' + '_empire' + '.txt'
        u_70Zn_empire2 = path + '/../EMPIRECOH2/' + foil + '/70Zn/' +  Z + '-' + 'Ni-66' + '_empire' + '.txt'
        u_70Zn_empire3 = path + '/../EMPIRECOH2/' + foil + '/70Zn/' +  Z + '-' + 'Cu-66' + '_empire' + '.txt'
        u_70Zn_empire4 = path + '/../EMPIRECOH2/' + foil + '/70Zn/' +  Z + '-' + 'Cu-67' + '_empire' + '.txt'
        u_70Zn_empire5 = path + '/../EMPIRECOH2/' + foil + '/70Zn/' +  Z + '-' + 'Zn-69' + '_empire' + '.txt'

        #f_66Zn = path + '/../EMPIRECOH2/' + foil + '/rp030066_' + Z +A + '_empire' + '.txt'
        #f_67Zn = path + '/../EMPIRECOH2/' + foil + '/rp030067_' + Z +A + file_ending + '.txt'
        #f_68Zn = path + '/../EMPIRECOH2/' + foil + '/rp030068_' + Z +A + file_ending + '.txt'
        #f_70Zn = path + '/../EMPIRECOH2/' + foil + '/rp030070_' + Z +A + file_ending + '.txt'


        print("Zn 64 file: ",u_64Zn_empire1)
        #print("Ir 191 file: ",f_191Ir)

############################ 64Zn #################################

        if os.path.isfile(u_64Zn_empire1): # Cu-64
            #print("Ir 191 file: ",f_191Ir)
            print("v_64Zn_Cu-64 exists")
            CS_64Zn_1 = np.genfromtxt(u_64Zn_empire1, delimiter='	', usecols=[1])
            E_64Zn_1 = np.genfromtxt(u_64Zn_empire1, delimiter='	', usecols=[0])
        else:
            print("64Zn_Cu-64 file does not exist")
            CS_64Zn_1 = 0
            E_64Zn_1 =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        #print('!!!!!!!!!!!!!!!!!!!!!!!!!')
        #print(E_64Zn_E)
        #plt.plot(E_64Zn_E, CS_64Zn_E, color='m', label='EMPIRE')

        if os.path.isfile(u_64Zn_empire2): # Zn-62
            #print("Ir 191 file: ",f_191Ir)
            print("v_64Zn_Zn-62 exists")
            CS_64Zn_2 = np.genfromtxt(u_64Zn_empire2, delimiter='	', usecols=[1])
            E_64Zn_2 = np.genfromtxt(u_64Zn_empire2, delimiter='	', usecols=[0])
        else:
            print("64Zn_Zn-62 file does not exist")
            CS_64Zn_2 = 0
            E_64Zn_2 =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        if os.path.isfile(u_64Zn_empire3): # Zn-63
            #print("Ir 191 file: ",f_191Ir)
            print("v_64Zn_Zn-63 exists")
            CS_64Zn_3 = np.genfromtxt(u_64Zn_empire3, delimiter='	', usecols=[1])
            E_64Zn_3 = np.genfromtxt(u_64Zn_empire3, delimiter='	', usecols=[0])
        else:
            print("64Zn_Zn-63 file does not exist")
            CS_64Zn_3 = 0
            E_64Zn_3 =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

############################ 66Zn #################################


        if os.path.isfile(u_66Zn_empire1): # Cu-64
            #print("Ir 191 file: ",f_191Ir)
            print("v_66Zn_Cu-64 exists")
            CS_66Zn_1 = np.genfromtxt(u_66Zn_empire1, delimiter='	', usecols=[1])
            E_66Zn_1 = np.genfromtxt(u_66Zn_empire1, delimiter='	', usecols=[0])
        else:
            print("66Zn_Cu-64 file does not exist")
            CS_66Zn_1 = 0
            E_66Zn_1 =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        if os.path.isfile(u_66Zn_empire2): # Cu-66
            #print("Ir 191 file: ",f_191Ir)
            print("v_66Zn_Cu-66 exists")
            CS_66Zn_2 = np.genfromtxt(u_66Zn_empire2, delimiter='	', usecols=[1])
            E_66Zn_2 = np.genfromtxt(u_66Zn_empire2, delimiter='	', usecols=[0])
        else:
            print("66Zn_Cu-66 file does not exist")
            CS_66Zn_2 = 0
            E_66Zn_2 =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        if os.path.isfile(u_66Zn_empire3): # Zn-65
            #print("Ir 191 file: ",f_191Ir)
            print("v_66Zn_Zn-65 exists")
            CS_66Zn_3 = np.genfromtxt(u_66Zn_empire3, delimiter='	', usecols=[1])
            E_66Zn_3 = np.genfromtxt(u_66Zn_empire3, delimiter='	', usecols=[0])
        else:
            print("66Zn_Zn-65 file does not exist")
            CS_66Zn_3 = 0
            E_66Zn_3 =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

############################ 67Zn #################################

        if os.path.isfile(u_67Zn_empire1): # Cu-67
            #print("Ir 191 file: ",f_191Ir)
            print("u_67Zn_Cu-67 exists")
            CS_67Zn_1 = np.genfromtxt(u_67Zn_empire1, delimiter='	', usecols=[1])
            E_67Zn_1 = np.genfromtxt(u_67Zn_empire1, delimiter='	', usecols=[0])
        else:
            print("67Zn_Cu-67 file does not exist")
            CS_67Zn_1 = 0
            E_67Zn_1 =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        #print('!!!!!!!!!!!!!!!!!!!!!!!!!')
        #print(E_68Zn_E)
        #plt.plot(E_68Zn_E, CS_68Zn_E, color='m', label='EMPIRE')


        if os.path.isfile(u_67Zn_empire2): # Ni-65
            #print("Ir 191 file: ",f_191Ir)
            print("u_67Zn_Ni-65 exists")
            CS_67Zn_2 = np.genfromtxt(u_67Zn_empire2, delimiter='	', usecols=[1])
            E_67Zn_2 = np.genfromtxt(u_67Zn_empire2, delimiter='	', usecols=[0])
        else:
            print("67Zn_Ni-65 file does not exist")
            CS_67Zn_2 = 0
            E_67Zn_2 =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        if os.path.isfile(u_67Zn_empire3): # Ni-66
            #print("Ir 191 file: ",f_191Ir)
            print("u_67Zn_Ni-66 exists")
            CS_67Zn_3 = np.genfromtxt(u_67Zn_empire3, delimiter='	', usecols=[1])
            E_67Zn_3 = np.genfromtxt(u_67Zn_empire3, delimiter='	', usecols=[0])
        else:
            print("67Zn_Ni-66 file does not exist")
            CS_67Zn_3 = 0
            E_67Zn_3 =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        if os.path.isfile(u_67Zn_empire4): # Cu-64
            #print("Ir 191 file: ",f_191Ir)
            print("u_67Zn_Cu-64 exists")
            CS_67Zn_4 = np.genfromtxt(u_67Zn_empire4, delimiter='	', usecols=[1])
            E_67Zn_4 = np.genfromtxt(u_67Zn_empire4, delimiter='	', usecols=[0])
        else:
            print("67Zn_Cu-64 file does not exist")
            CS_67Zn_4 = 0
            E_67Zn_4 =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        if os.path.isfile(u_67Zn_empire5): # Cu-66
            #print("Ir 191 file: ",f_191Ir)
            print("u_67Zn_Cu-66 exists")
            CS_67Zn_5 = np.genfromtxt(u_67Zn_empire5, delimiter='	', usecols=[1])
            E_67Zn_5 = np.genfromtxt(u_67Zn_empire5, delimiter='	', usecols=[0])
        else:
            print("67Zn_Cu-66 file does not exist")
            CS_67Zn_5 = 0
            E_67Zn_5 =  0


        if os.path.isfile(u_67Zn_empire6): # Cu-67
            #print("Ir 191 file: ",f_191Ir)
            print("u_67Zn_Cu-67 exists")
            CS_67Zn_6 = np.genfromtxt(u_67Zn_empire6, delimiter='	', usecols=[1])
            E_67Zn_6 = np.genfromtxt(u_67Zn_empire6, delimiter='	', usecols=[0])
        else:
            print("67Zn_Cu-67 file does not exist")
            CS_67Zn_6 = 0
            E_67Zn_6 =  0


        if os.path.isfile(u_67Zn_empire7): # Zn-65
            #print("Ir 191 file: ",f_191Ir)
            print("u_67Zn_Zn-65 exists")
            CS_67Zn_7 = np.genfromtxt(u_67Zn_empire7, delimiter='	', usecols=[1])
            E_67Zn_7 = np.genfromtxt(u_67Zn_empire7, delimiter='	', usecols=[0])
        else:
            print("67Zn_Zn-65 file does not exist")
            CS_67Zn_7 = 0
            E_67Zn_7 =  0

############################ 68Zn #################################

        if os.path.isfile(u_68Zn_empire1): # Ni-65
            #print("Ir 191 file: ",f_191Ir)
            print("u_68Zn_Ni-65 exists")
            CS_68Zn_1 = np.genfromtxt(u_68Zn_empire1, delimiter='	', usecols=[1])
            E_68Zn_1 = np.genfromtxt(u_68Zn_empire1, delimiter='	', usecols=[0])
        else:
            print("68Zn_Ni-65 file does not exist")
            CS_68Zn_1 = 0
            E_68Zn_1 =  0


        if os.path.isfile(u_68Zn_empire2): # Ni-66
            #print("Ir 191 file: ",f_191Ir)
            print("u_68Zn_Ni-66 exists")
            CS_68Zn_2 = np.genfromtxt(u_68Zn_empire2, delimiter='	', usecols=[1])
            E_68Zn_2 = np.genfromtxt(u_68Zn_empire2, delimiter='	', usecols=[0])
        else:
            print("68Zn_Ni-66 file does not exist")
            CS_68Zn_2 = 0
            E_68Zn_2 =  0


        if os.path.isfile(u_68Zn_empire3): # Cu-64
            #print("Ir 191 file: ",f_191Ir)
            print("u_68Zn_Cu-64 exists")
            CS_68Zn_3 = np.genfromtxt(u_68Zn_empire3, delimiter='	', usecols=[1])
            E_68Zn_3 = np.genfromtxt(u_68Zn_empire3, delimiter='	', usecols=[0])
        else:
            print("68Zn_Cu-64 file does not exist")
            CS_68Zn_3 = 0
            E_68Zn_3 =  0


        if os.path.isfile(u_68Zn_empire4): # Cu-66
            #print("Ir 191 file: ",f_191Ir)
            print("u_68Zn_Cu-64 exists")
            CS_68Zn_4 = np.genfromtxt(u_68Zn_empire4, delimiter='	', usecols=[1])
            E_68Zn_4 = np.genfromtxt(u_68Zn_empire4, delimiter='	', usecols=[0])
        else:
            print("68Zn_Cu-64 file does not exist")
            CS_68Zn_4 = 0
            E_68Zn_4 =  0


        if os.path.isfile(u_68Zn_empire5): # Cu-67
            #print("Ir 191 file: ",f_191Ir)
            print("u_68Zn_Cu-67 exists")
            CS_68Zn_5 = np.genfromtxt(u_68Zn_empire5, delimiter='	', usecols=[1])
            E_68Zn_5 = np.genfromtxt(u_68Zn_empire5, delimiter='	', usecols=[0])
        else:
            print("68Zn_Cu-67 file does not exist")
            CS_68Zn_5 = 0
            E_68Zn_5 =  0


        if os.path.isfile(u_68Zn_empire6): # Zn-65
            #print("Ir 191 file: ",f_191Ir)
            print("u_68Zn_Zn-65 exists")
            CS_68Zn_6 = np.genfromtxt(u_68Zn_empire6, delimiter='	', usecols=[1])
            E_68Zn_6 = np.genfromtxt(u_68Zn_empire6, delimiter='	', usecols=[0])
        else:
            print("68Zn_Zn-65 file does not exist")
            CS_68Zn_6 = 0
            E_68Zn_6 =  0

############################ 70Zn #################################

        if os.path.isfile(u_70Zn_empire1): # Ni-65
            #print("Ir 191 file: ",f_191Ir)
            print("u_70Zn_Ni-65 exists")
            CS_70Zn_1 = np.genfromtxt(u_70Zn_empire1, delimiter='	', usecols=[1])
            E_70Zn_1 = np.genfromtxt(u_70Zn_empire1, delimiter='	', usecols=[0])
        else:
            print("70Zn_Ni-65 file does not exist")
            CS_70Zn_1 = 0
            E_70Zn_1 =  0


        if os.path.isfile(u_70Zn_empire2): # Ni-66
            #print("Ir 191 file: ",f_191Ir)
            print("u_70Zn_Ni-66 exists")
            CS_70Zn_2 = np.genfromtxt(u_70Zn_empire2, delimiter='	', usecols=[1])
            E_70Zn_2 = np.genfromtxt(u_70Zn_empire2, delimiter='	', usecols=[0])
        else:
            print("70Zn_Ni-66 file does not exist")
            CS_70Zn_2 = 0
            E_70Zn_2 =  0


        if os.path.isfile(u_70Zn_empire3): # Cu-66
            #print("Ir 191 file: ",f_191Ir)
            print("u_70Zn_Cu-66 exists")
            CS_70Zn_3 = np.genfromtxt(u_70Zn_empire3, delimiter='	', usecols=[1])
            E_70Zn_3 = np.genfromtxt(u_70Zn_empire3, delimiter='	', usecols=[0])
        else:
            print("70Zn_Cu-66 file does not exist")
            CS_70Zn_3 = 0
            E_70Zn_3 =  0


        if os.path.isfile(u_70Zn_empire4): # Cu-67
            #print("Ir 191 file: ",f_191Ir)
            print("u_70Zn_Cu-67 exists")
            CS_70Zn_4 = np.genfromtxt(u_70Zn_empire4, delimiter='	', usecols=[1])
            E_70Zn_4 = np.genfromtxt(u_70Zn_empire4, delimiter='	', usecols=[0])
        else:
            print("70Zn_Cu-67 file does not exist")
            CS_70Zn_4 = 0
            E_70Zn_4 =  0


        if os.path.isfile(u_70Zn_empire5): # Zn-69
            #print("Ir 191 file: ",f_191Ir)
            print("u_70Zn_Zn-69 exists")
            CS_70Zn_5 = np.genfromtxt(u_70Zn_empire5, delimiter='	', usecols=[1])
            E_70Zn_5 = np.genfromtxt(u_70Zn_empire5, delimiter='	', usecols=[0])
        else:
            print("70Zn_Zn-69 file does not exist")
            CS_70Zn_5 = 0
            E_70Zn_5 =  0




        E_64Zn = E_64Zn_1*abund_64Zn + E_64Zn_2*abund_64Zn + E_64Zn_3*abund_64Zn
        E_66Zn = E_66Zn_1*abund_66Zn + E_66Zn_2*abund_66Zn + E_66Zn_3*abund_66Zn
        E_67Zn = E_67Zn_1*abund_67Zn + E_67Zn_2*abund_67Zn + E_67Zn_3*abund_67Zn + E_67Zn_4*abund_67Zn + E_67Zn_5*abund_67Zn + E_67Zn_6*abund_67Zn + E_67Zn_7*abund_67Zn
        E_68Zn = E_68Zn_1*abund_68Zn + E_68Zn_2*abund_68Zn + E_68Zn_3*abund_68Zn + E_68Zn_4*abund_68Zn + E_68Zn_5*abund_68Zn + E_68Zn_6*abund_68Zn
        E_70Zn = E_70Zn_1*abund_70Zn + E_70Zn_2*abund_70Zn + E_70Zn_3*abund_70Zn + E_70Zn_4*abund_70Zn + E_70Zn_5*abund_70Zn

        E = E_64Zn + E_66Zn + E_67Zn + E_68Zn + E_70Zn


        CS_64Zn = CS_64Zn_1*abund_64Zn + CS_64Zn_2*abund_64Zn + CS_64Zn_3*abund_64Zn
        CS_66Zn = CS_66Zn_1*abund_66Zn + CS_66Zn_2*abund_66Zn + CS_66Zn_3*abund_66Zn
        CS_67Zn = CS_67Zn_1*abund_67Zn + CS_67Zn_2*abund_67Zn + CS_67Zn_3*abund_67Zn + CS_67Zn_4*abund_67Zn + CS_67Zn_5*abund_67Zn + CS_67Zn_6*abund_67Zn + CS_67Zn_7*abund_67Zn
        CS_68Zn = CS_68Zn_1*abund_68Zn + CS_68Zn_2*abund_68Zn + CS_68Zn_3*abund_68Zn + CS_68Zn_4*abund_68Zn + CS_68Zn_5*abund_68Zn + CS_68Zn_6*abund_68Zn
        CS_70Zn = CS_70Zn_1*abund_70Zn + CS_70Zn_2*abund_70Zn + CS_70Zn_3*abund_70Zn + CS_70Zn_4*abund_70Zn + CS_70Zn_5*abund_70Zn


        E = E_64Zn + E_66Zn + E_67Zn + E_68Zn + E_70Zn
        CS = CS_64Zn + CS_66Zn + CS_67Zn + CS_68Zn + CS_70Zn

        #
        #print(len(E))
        #print(len(CS))
        print('...........................')
        #plt.plot(E, xaxis, 'y--', label='TENDL')
        plt.plot(E, CS, 'y--', label='EMPIRE')

#--------------------- Al ------------------------

    if foil == 'Al':
    #A = ['191', '193'] # stable iridium isotopes
        abund_27Al = 1.00
    #file_ending =
    #endre f_191I til de stabile isotopene i zink osv
    #v_64Zn = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-0' + filename_CS + '_coh' + '.txt'
    # 64Zn
        u_27Al_empire = path + '/../EMPIRECOH2/' + foil + '/27Al/' +  Z + '-' + 'Na-24' + '_empire' + '.txt'


        if os.path.isfile(u_27Al_empire): # Na-24
            #print("Ir 191 file: ",f_191Ir)
            print("u_27Al_Na-24 exists")
            CS_27Al = np.genfromtxt(u_27Al_empire, delimiter='	', usecols=[1])
            E_27Al = np.genfromtxt(u_27Al_empire, delimiter='	', usecols=[0])
        else:
            print("27Al_Na-24 file does not exist")
            CS_27Al = 0
            E_27Al =  0

        E = E_27Al*abund_27Al
        CS = CS_27Al*abund_27Al
        plt.plot(E, CS, 'y--', label='EMPIRE')

#------------------------- In -----------------------------

    if foil == 'In':
    #A = ['191', '193'] # stable iridium isotopes
        abund_113In = 0.429 ; abund_115In = 0.9571
    #file_ending =
    #endre f_191I til de stabile isotopene i zink osv
    #v_64Zn = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-0' + filename_CS + '_coh' + '.txt'
        # 113In
        u_113In_empire1 = path + '/../EMPIRECOH2/' + foil + '/113In/' +  Z + '-' + 'In-111' + '_empire' + '.txt'
        u_113In_empire2 = path + '/../EMPIRECOH2/' + foil + '/113In/' +  Z + '-' + 'In-112' + '_empire' + '.txt'
        u_113In_empire3 = path + '/../EMPIRECOH2/' + foil + '/113In/' +  Z + '-' + 'In-112M' + '_empire' + '.txt'
        # 115In
        u_115In_empire1 = path + '/../EMPIRECOH2/' + foil + '/115In/' +  Z + '-' + 'In-112' + '_empire' + '.txt'
        u_115In_empire2 = path + '/../EMPIRECOH2/' + foil + '/115In/' +  Z + '-' + 'In-112M' + '_empire' + '.txt'
        u_115In_empire3 = path + '/../EMPIRECOH2/' + foil + '/115In/' +  Z + '-' + 'In-113' + '_empire' + '.txt'
        u_115In_empire4 = path + '/../EMPIRECOH2/' + foil + '/115In/' +  Z + '-' + 'In-113M' + '_empire' + '.txt'
        u_115In_empire5 = path + '/../EMPIRECOH2/' + foil + '/115In/' +  Z + '-' + 'In-114' + '_empire' + '.txt'
        u_115In_empire6 = path + '/../EMPIRECOH2/' + foil + '/115In/' +  Z + '-' + 'In-114M' + '_empire' + '.txt'


########################## 113In #############################

        if os.path.isfile(u_113In_empire1): # In-111
            #print("Ir 191 file: ",f_191Ir)
            print("u_113In_In-111 exists")
            CS_113In_1 = np.genfromtxt(u_113In_empire1, delimiter='	', usecols=[1])
            E_113In_1 = np.genfromtxt(u_113In_empire1, delimiter='	', usecols=[0])
        else:
            print("113In_In-111 file does not exist")
            CS_113In_1 = 0
            E_113In_1 =  0


        if os.path.isfile(u_113In_empire2): # In-112
            #print("Ir 191 file: ",f_191Ir)
            print("u_113In_In-112 exists")
            CS_113In_2 = np.genfromtxt(u_113In_empire2, delimiter='	', usecols=[1])
            E_113In_2 = np.genfromtxt(u_113In_empire2, delimiter='	', usecols=[0])
        else:
            print("113In_In-112 file does not exist")
            CS_113In_2 = 0
            E_113In_2 =  0


        if os.path.isfile(u_113In_empire3): # In-112M
            #print("Ir 191 file: ",f_191Ir)
            print("u_113In_In-112M exists")
            CS_113In_3 = np.genfromtxt(u_113In_empire3, delimiter='	', usecols=[1])
            E_113In_3 = np.genfromtxt(u_113In_empire3, delimiter='	', usecols=[0])
        else:
            print("113In_In-112M file does not exist")
            CS_113In_3 = 0
            E_113In_3 =  0


########################## 115In #############################

        if os.path.isfile(u_115In_empire1): # In-112
            #print("Ir 191 file: ",f_191Ir)
            print("u_115In_In-112 exists")
            CS_115In_1 = np.genfromtxt(u_115In_empire1, delimiter='	', usecols=[1])
            E_115In_1 = np.genfromtxt(u_115In_empire1, delimiter='	', usecols=[0])
        else:
            print("115In_In-112 file does not exist")
            CS_115In_1 = 0
            E_115In_1 =  0


        if os.path.isfile(u_115In_empire2): # In-112M
            #print("Ir 191 file: ",f_191Ir)
            print("u_115In_In-112M exists")
            CS_115In_2 = np.genfromtxt(u_115In_empire2, delimiter='	', usecols=[1])
            E_115In_2 = np.genfromtxt(u_115In_empire2, delimiter='	', usecols=[0])
        else:
            print("115In_In-112M file does not exist")
            CS_115In_2 = 0
            E_115In_2 =  0


        if os.path.isfile(u_115In_empire3): # In-113
            #print("Ir 191 file: ",f_191Ir)
            print("u_115In_In-113 exists")
            CS_115In_3 = np.genfromtxt(u_115In_empire3, delimiter='	', usecols=[1])
            E_115In_3 = np.genfromtxt(u_115In_empire3, delimiter='	', usecols=[0])
        else:
            print("115In_In-113 file does not exist")
            CS_115In_3 = 0
            E_113In_3 =  0


        if os.path.isfile(u_115In_empire4): # In-113M
            #print("Ir 191 file: ",f_191Ir)
            print("u_115In_In-113M exists")
            CS_115In_4 = np.genfromtxt(u_115In_empire4, delimiter='	', usecols=[1])
            E_115In_4 = np.genfromtxt(u_115In_empire4, delimiter='	', usecols=[0])
        else:
            print("115In_In-113M file does not exist")
            CS_115In_4 = 0
            E_115In_4 =  0


        if os.path.isfile(u_115In_empire5): # In-114
            #print("Ir 191 file: ",f_191Ir)
            print("u_115In_In-114 exists")
            CS_115In_5 = np.genfromtxt(u_115In_empire5, delimiter='	', usecols=[1])
            E_115In_5 = np.genfromtxt(u_115In_empire5, delimiter='	', usecols=[0])
        else:
            print("115In_In-114 file does not exist")
            CS_115In_5 = 0
            E_115In_5 =  0


        if os.path.isfile(u_115In_empire6): # In-114M
            #print("Ir 191 file: ",f_191Ir)
            print("u_115In_In-114M exists")
            CS_115In_6 = np.genfromtxt(u_115In_empire6, delimiter='	', usecols=[1])
            E_115In_6 = np.genfromtxt(u_115In_empire6, delimiter='	', usecols=[0])
        else:
            print("115In_In-114M file does not exist")
            CS_115In_6 = 0
            E_115In_6 =  0

        E_113In = E_113In_1*abund_113In + E_113In_2*abund_113In
        E_115In = E_115In_1*abund_115In + E_115In_2*abund_115In + E_115In_3*abund_115In + E_115In_4*abund_115In+ E_115In_5*abund_115In+ E_115In_6*abund_115In

        CS_113In = (CS_113In_1 + CS_113In_2 ) * abund_113In
        CS_115In = (CS_115In_1 + CS_115In_2 + CS_115In_3 + CS_115In_4 + CS_115In_5 + CS_115In_6) * abund_115In

        E = E_113In +  E_115In
        CS = CS_113In + CS_115In

        plt.plot(E, CS, 'y--', label='EMPIRE')


#------------------------- Y -----------------------------

    if foil == 'Y':
    #A = ['191', '193'] # stable iridium isotopes
        abund_89Y = 1.00 ;
    #file_ending =
    #endre f_191I til de stabile isotopene i zink osv
    #v_64Zn = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-0' + filename_CS + '_coh' + '.txt'
        # 89Y
        u_89Y_empire1 = path + '/../EMPIRECOH2/' + foil + '/89Y/' +  Z + '-' + 'Y-87' + '_empire' + '.txt'
        u_89Y_empire2 = path + '/../EMPIRECOH2/' + foil + '/89Y/' +  Z + '-' + 'Y-88' + '_empire' + '.txt'


        if os.path.isfile(u_89Y_empire1): # Y-87
            #print("Ir 191 file: ",f_191Ir)
            print("u_89Y_Y-87 exists")
            CS_89Y_1 = np.genfromtxt(u_89Y_empire1, delimiter='	', usecols=[1])
            E_89Y_1 = np.genfromtxt(u_89Y_empire1, delimiter='	', usecols=[0])
        else:
            print("89Y_Y-87 file does not exist")
            CS_89Y_1 = 0
            E_89Y_1 =  0


        if os.path.isfile(u_89Y_empire2): # Y-88
            #print("Ir 191 file: ",f_191Ir)
            print("u_89Y_Y-88 exists")
            CS_89Y_2 = np.genfromtxt(u_89Y_empire2, delimiter='	', usecols=[1])
            E_89Y_2 = np.genfromtxt(u_89Y_empire2, delimiter='	', usecols=[0])
        else:
            print("89Y_Y-88 file does not exist")
            CS_89Y_2 = 0
            E_89Y_2 =  0


        E = (E_89Y_1 + E_89Y_2) *  abund_89Y
        CS = (CS_89Y_1 + CS_89Y_2) * abund_89Y
        plt.plot(E, CS, 'y--', label='EMPIRE')

#------------------------- Zr -----------------------------


    if foil == 'Zr':
    #A = ['191', '193'] # stable iridium isotopes
        abund_90Zr = 0.5145 ; abund_91Zr = 0.1122 ; abund_92Zr = 0.1715 ; abund_94Zr = 0.1738 ; abund_96Zr = 0.280
    #file_ending =
    #endre f_191I til de stabile isotopene i zink osv
    #v_64Zn = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-0' + filename_CS + '_coh' + '.txt'

######################### 90 Zr ############################
        u_90Zr_empire1 = path + '/../EMPIRECOH2/' + foil + '/90Zr/' +  Z + '-' + 'Y-90' + '_empire' + '.txt'
        u_90Zr_empire2 = path + '/../EMPIRECOH2/' + foil + '/90Zr/' +  Z + '-' + 'Y-90M' + '_empire' + '.txt'
        u_90Zr_empire3 = path + '/../EMPIRECOH2/' + foil + '/90Zr/' +  Z + '-' + 'Zr-89' + '_empire' + '.txt'

        if os.path.isfile(u_90Zr_empire1): # Y-90
            #print("Ir 191 file: ",f_191Ir)
            print("u_90Zr_Y-90 exists")
            CS_90Zr_1 = np.genfromtxt(u_90Zr_empire1, delimiter='	', usecols=[1])
            E_90Zr_1 = np.genfromtxt(u_90Zr_empire1, delimiter='	', usecols=[0])
        else:
            print("90Zr_Y-90 file does not exist")
            CS_90Zr_1 = 0
            E_90Zr_1 =  0

        if os.path.isfile(u_90Zr_empire2): # Y-90m
            #print("Ir 191 file: ",f_191Ir)
            print("u_90Zr_Y-90m exists")
            CS_90Zr_2 = np.genfromtxt(u_90Zr_empire2, delimiter='	', usecols=[1])
            E_90Zr_2 = np.genfromtxt(u_90Zr_empire2, delimiter='	', usecols=[0])
        else:
            print("90Zr_Y-90m file does not exist")
            CS_90Zr_2 = 0
            E_90Zr_2 =  0

        if os.path.isfile(u_90Zr_empire3): # Zr-89
            #print("Ir 191 file: ",f_191Ir)
            print("u_90Zr_Zr-89 exists")
            CS_90Zr_3 = np.genfromtxt(u_90Zr_empire3, delimiter='	', usecols=[1])
            E_90Zr_3 = np.genfromtxt(u_90Zr_empire3, delimiter='	', usecols=[0])
        else:
            print("90Zr_Zr-89 file does not exist")
            CS_90Zr_3 = 0
            E_90Zr_3 =  0


######################### 91 Zr ############################
        u_91Zr_empire1 = path + '/../EMPIRECOH2/' + foil + '/91Zr/' +  Z + '-' + 'Y-90' + '_empire' + '.txt'
        u_91Zr_empire2 = path + '/../EMPIRECOH2/' + foil + '/91Zr/' +  Z + '-' + 'Y-90M' + '_empire' + '.txt'
        u_91Zr_empire3 = path + '/../EMPIRECOH2/' + foil + '/91Zr/' +  Z + '-' + 'Y-91' + '_empire' + '.txt'
        u_91Zr_empire4 = path + '/../EMPIRECOH2/' + foil + '/91Zr/' +  Z + '-' + 'Y-91M' + '_empire' + '.txt'
        u_91Zr_empire5 = path + '/../EMPIRECOH2/' + foil + '/91Zr/' +  Z + '-' + 'Zr-89' + '_empire' + '.txt'


        if os.path.isfile(u_91Zr_empire1): # Y-90
            #print("Ir 191 file: ",f_191Ir)
            print("u_91Zr_Y-90 exists")
            CS_91Zr_1 = np.genfromtxt(u_91Zr_empire1, delimiter='	', usecols=[1])
            E_91Zr_1 = np.genfromtxt(u_91Zr_empire1, delimiter='	', usecols=[0])
        else:
            print("91Zr_Y-90 file does not exist")
            CS_91Zr_1 = 0
            E_91Zr_1 =  0

        if os.path.isfile(u_91Zr_empire2): # Y-90M
            #print("Ir 191 file: ",f_191Ir)
            print("u_91Zr_Y-90M exists")
            CS_91Zr_2 = np.genfromtxt(u_91Zr_empire2, delimiter='	', usecols=[1])
            E_91Zr_2 = np.genfromtxt(u_91Zr_empire2, delimiter='	', usecols=[0])
        else:
            print("91Zr_Y-90M file does not exist")
            CS_91Zr_2 = 0
            E_91Zr_2 =  0

        if os.path.isfile(u_91Zr_empire3): # Y-91
            #print("Ir 191 file: ",f_191Ir)
            print("u_91Zr_Y-91 exists")
            CS_91Zr_3 = np.genfromtxt(u_91Zr_empire3, delimiter='	', usecols=[1])
            E_91Zr_3 = np.genfromtxt(u_91Zr_empire3, delimiter='	', usecols=[0])
        else:
            print("91Zr_Y-91 file does not exist")
            CS_91Zr_3 = 0
            E_91Zr_3 =  0

        if os.path.isfile(u_91Zr_empire4): # Y-91M
            #print("Ir 191 file: ",f_191Ir)
            print("u_91Zr_Y-91M exists")
            CS_91Zr_4 = np.genfromtxt(u_91Zr_empire4, delimiter='	', usecols=[1])
            E_91Zr_4 = np.genfromtxt(u_91Zr_empire4, delimiter='	', usecols=[0])
        else:
            print("91Zr_Y-91M file does not exist")
            CS_91Zr_4 = 0
            E_91Zr_4 =  0

        if os.path.isfile(u_91Zr_empire5): # Zr-89
            #print("Ir 191 file: ",f_191Ir)
            print("u_91Zr_Zr-89 exists")
            CS_91Zr_5 = np.genfromtxt(u_91Zr_empire5, delimiter='	', usecols=[1])
            E_91Zr_5 = np.genfromtxt(u_91Zr_empire5, delimiter='	', usecols=[0])
        else:
            print("91Zr_Zr-89 file does not exist")
            CS_91Zr_5 = 0
            E_91Zr_5 =  0


######################### 92 Zr ############################
        u_92Zr_empire1 = path + '/../EMPIRECOH2/' + foil + '/92Zr/' +  Z + '-' + 'Y-90' + '_empire' + '.txt'
        u_92Zr_empire2 = path + '/../EMPIRECOH2/' + foil + '/92Zr/' +  Z + '-' + 'Y-90M' + '_empire' + '.txt'
        u_92Zr_empire3 = path + '/../EMPIRECOH2/' + foil + '/92Zr/' +  Z + '-' + 'Y-91' + '_empire' + '.txt'
        u_92Zr_empire4 = path + '/../EMPIRECOH2/' + foil + '/92Zr/' +  Z + '-' + 'Y-91M' + '_empire' + '.txt'
        u_92Zr_empire5 = path + '/../EMPIRECOH2/' + foil + '/92Zr/' +  Z + '-' + 'Y-92' + '_empire' + '.txt'
        u_92Zr_empire6 = path + '/../EMPIRECOH2/' + foil + '/92Zr/' +  Z + '-' + 'Zr-89' + '_empire' + '.txt'

        if os.path.isfile(u_92Zr_empire1): # Y-90
            #print("Ir 191 file: ",f_191Ir)
            print("u_92Zr_Y-90 exists")
            CS_92Zr_1 = np.genfromtxt(u_92Zr_empire1, delimiter='	', usecols=[1])
            E_92Zr_1 = np.genfromtxt(u_92Zr_empire1, delimiter='	', usecols=[0])
        else:
            print("92Zr_Y-90 file does not exist")
            CS_92Zr_1 = 0
            E_92Zr_1 =  0

        if os.path.isfile(u_92Zr_empire2): # Y-90M
            #print("Ir 191 file: ",f_191Ir)
            print("u_92Zr_Y-90M exists")
            CS_92Zr_2 = np.genfromtxt(u_92Zr_empire2, delimiter='	', usecols=[1])
            E_92Zr_2 = np.genfromtxt(u_92Zr_empire2, delimiter='	', usecols=[0])
        else:
            print("92Zr_Y-90M file does not exist")
            CS_92Zr_2 = 0
            E_92Zr_2 =  0

        if os.path.isfile(u_92Zr_empire3): # Y-91
            #print("Ir 191 file: ",f_191Ir)
            print("u_92Zr_Y-91 exists")
            CS_92Zr_3 = np.genfromtxt(u_92Zr_empire3, delimiter='	', usecols=[1])
            E_92Zr_3 = np.genfromtxt(u_92Zr_empire3, delimiter='	', usecols=[0])
        else:
            print("92Zr_Y-91 file does not exist")
            CS_92Zr_3 = 0
            E_92Zr_3 =  0

        if os.path.isfile(u_92Zr_empire4): # Y-91M
            #print("Ir 191 file: ",f_191Ir)
            print("u_92Zr_Y-91M exists")
            CS_92Zr_4 = np.genfromtxt(u_92Zr_empire4, delimiter='	', usecols=[1])
            E_92Zr_4 = np.genfromtxt(u_92Zr_empire4, delimiter='	', usecols=[0])
        else:
            print("92Zr_Y-91M file does not exist")
            CS_92Zr_4 = 0
            E_92Zr_4 =  0

        if os.path.isfile(u_92Zr_empire5): # Y-92
            #print("Ir 191 file: ",f_191Ir)
            print("u_92Zr_Y-92 exists")
            CS_92Zr_5 = np.genfromtxt(u_92Zr_empire5, delimiter='	', usecols=[1])
            E_92Zr_5 = np.genfromtxt(u_92Zr_empire5, delimiter='	', usecols=[0])
        else:
            print("92Zr_Y-92 file does not exist")
            CS_92Zr_5 = 0
            E_92Zr_5 =  0

        if os.path.isfile(u_92Zr_empire6): # Zr-89
            #print("Ir 191 file: ",f_191Ir)
            print("u_92Zr_Zr-89 exists")
            CS_92Zr_6 = np.genfromtxt(u_92Zr_empire6, delimiter='	', usecols=[1])
            E_92Zr_6 = np.genfromtxt(u_92Zr_empire6, delimiter='	', usecols=[0])
        else:
            print("92Zr_Zr-89 file does not exist")
            CS_92Zr_6 = 0
            E_92Zr_6 =  0


######################### 94 Zr ############################
        u_94Zr_empire1 = path + '/../EMPIRECOH2/' + foil + '/94Zr/' +  Z + '-' + 'Sr-91' + '_empire' + '.txt'
        u_94Zr_empire2 = path + '/../EMPIRECOH2/' + foil + '/94Zr/' +  Z + '-' + 'Sr-92' + '_empire' + '.txt'
        u_94Zr_empire3 = path + '/../EMPIRECOH2/' + foil + '/94Zr/' +  Z + '-' + 'Y-90' + '_empire' + '.txt'
        u_94Zr_empire4 = path + '/../EMPIRECOH2/' + foil + '/94Zr/' +  Z + '-' + 'Y-90M' + '_empire' + '.txt'
        u_94Zr_empire5 = path + '/../EMPIRECOH2/' + foil + '/94Zr/' +  Z + '-' + 'Y-91' + '_empire' + '.txt'
        u_94Zr_empire6 = path + '/../EMPIRECOH2/' + foil + '/94Zr/' +  Z + '-' + 'Y-91M' + '_empire' + '.txt'
        u_94Zr_empire7 = path + '/../EMPIRECOH2/' + foil + '/94Zr/' +  Z + '-' + 'Y-92' + '_empire' + '.txt'

        if os.path.isfile(u_94Zr_empire1): # Sr-91
            #print("Ir 191 file: ",f_191Ir)
            print("u_94Zr_Sr-91 exists")
            CS_94Zr_1 = np.genfromtxt(u_94Zr_empire1, delimiter='	', usecols=[1])
            E_94Zr_1 = np.genfromtxt(u_94Zr_empire1, delimiter='	', usecols=[0])
        else:
            print("94Zr_Sr-91 file does not exist")
            CS_94Zr_1 = 0
            E_94Zr_1 =  0

        if os.path.isfile(u_94Zr_empire2): # Sr-92
            #print("Ir 191 file: ",f_191Ir)
            print("u_94Zr_Sr-92 exists")
            CS_94Zr_2 = np.genfromtxt(u_94Zr_empire2, delimiter='	', usecols=[1])
            E_94Zr_2 = np.genfromtxt(u_94Zr_empire2, delimiter='	', usecols=[0])
        else:
            print("94Zr_Sr-92 file does not exist")
            CS_94Zr_2 = 0
            E_94Zr_2 =  0

        if os.path.isfile(u_94Zr_empire3): # Y-90
            #print("Ir 191 file: ",f_191Ir)
            print("u_94Zr_Y-90 exists")
            CS_94Zr_3 = np.genfromtxt(u_94Zr_empire3, delimiter='	', usecols=[1])
            E_94Zr_3 = np.genfromtxt(u_94Zr_empire3, delimiter='	', usecols=[0])
        else:
            print("94Zr_Y-90 file does not exist")
            CS_94Zr_3 = 0
            E_94Zr_3 =  0

        if os.path.isfile(u_94Zr_empire4): # Y-90M
            #print("Ir 191 file: ",f_191Ir)
            print("u_94Zr_Y-90M exists")
            CS_94Zr_4 = np.genfromtxt(u_94Zr_empire4, delimiter='	', usecols=[1])
            E_94Zr_4 = np.genfromtxt(u_94Zr_empire4, delimiter='	', usecols=[0])
        else:
            print("94Zr_Y-90M file does not exist")
            CS_94Zr_4 = 0
            E_94Zr_4 =  0

        if os.path.isfile(u_94Zr_empire5): # Y-91
            #print("Ir 191 file: ",f_191Ir)
            print("u_94Zr_Y-91 exists")
            CS_94Zr_5 = np.genfromtxt(u_94Zr_empire5, delimiter='	', usecols=[1])
            E_94Zr_5 = np.genfromtxt(u_94Zr_empire5, delimiter='	', usecols=[0])
        else:
            print("94Zr_Y-91 file does not exist")
            CS_94Zr_5 = 0
            E_94Zr_5 =  0

        if os.path.isfile(u_94Zr_empire6): # Y-91M
            #print("Ir 191 file: ",f_191Ir)
            print("u_94Zr_Y-91M exists")
            CS_94Zr_6 = np.genfromtxt(u_94Zr_empire6, delimiter='	', usecols=[1])
            E_94Zr_6 = np.genfromtxt(u_94Zr_empire6, delimiter='	', usecols=[0])
        else:
            print("94Zr_Y-91M file does not exist")
            CS_94Zr_6 = 0
            E_94Zr_6 =  0

        if os.path.isfile(u_94Zr_empire7): # Y-92
            #print("Ir 191 file: ",f_191Ir)
            print("u_94Zr_Y-92 exists")
            CS_94Zr_7 = np.genfromtxt(u_94Zr_empire7, delimiter='	', usecols=[1])
            E_94Zr_7 = np.genfromtxt(u_94Zr_empire7, delimiter='	', usecols=[0])
        else:
            print("94Zr_Y-92 file does not exist")
            CS_94Zr_7 = 0
            E_94Zr_7 =  0


######################### 96 Zr ############################
        u_96Zr_empire1 = path + '/../EMPIRECOH2/' + foil + '/96Zr/' +  Z + '-' + 'Sr-91' + '_empire' + '.txt'
        u_96Zr_empire2 = path + '/../EMPIRECOH2/' + foil + '/96Zr/' +  Z + '-' + 'Sr-92' + '_empire' + '.txt'
        u_96Zr_empire3 = path + '/../EMPIRECOH2/' + foil + '/96Zr/' +  Z + '-' + 'Sr-94' + '_empire' + '.txt'
        u_96Zr_empire4 = path + '/../EMPIRECOH2/' + foil + '/96Zr/' +  Z + '-' + 'Y-91' + '_empire' + '.txt'
        u_96Zr_empire5 = path + '/../EMPIRECOH2/' + foil + '/96Zr/' +  Z + '-' + 'Y-91M' + '_empire' + '.txt'
        u_96Zr_empire6 = path + '/../EMPIRECOH2/' + foil + '/96Zr/' +  Z + '-' + 'Y-92' + '_empire' + '.txt'
        u_96Zr_empire7 = path + '/../EMPIRECOH2/' + foil + '/96Zr/' +  Z + '-' + 'Y-95' + '_empire' + '.txt'
        u_96Zr_empire8 = path + '/../EMPIRECOH2/' + foil + '/96Zr/' +  Z + '-' + 'Zr-95' + '_empire' + '.txt'

        if os.path.isfile(u_96Zr_empire1): # Sr-91
            #print("Ir 191 file: ",f_191Ir)
            print("u_96Zr_Sr-91 exists")
            CS_96Zr_1 = np.genfromtxt(u_96Zr_empire1, delimiter='	', usecols=[1])
            E_96Zr_1 = np.genfromtxt(u_96Zr_empire1, delimiter='	', usecols=[0])
        else:
            print("96Zr_Sr-91 file does not exist")
            CS_96Zr_1 = 0
            E_96Zr_1 =  0

        if os.path.isfile(u_96Zr_empire2): # Sr-92
            #print("Ir 191 file: ",f_191Ir)
            print("u_96Zr_Sr-92 exists")
            CS_96Zr_2 = np.genfromtxt(u_96Zr_empire2, delimiter='	', usecols=[1])
            E_96Zr_2 = np.genfromtxt(u_96Zr_empire2, delimiter='	', usecols=[0])
        else:
            print("96Zr_Sr-92 file does not exist")
            CS_96Zr_2 = 0
            E_96Zr_2 =  0

        if os.path.isfile(u_96Zr_empire3): # Sr-94
            #print("Ir 191 file: ",f_191Ir)
            print("u_96Zr_Sr-94 exists")
            CS_96Zr_3 = np.genfromtxt(u_96Zr_empire3, delimiter='	', usecols=[1])
            E_96Zr_3 = np.genfromtxt(u_96Zr_empire3, delimiter='	', usecols=[0])
        else:
            print("96Zr_Sr-94 file does not exist")
            CS_96Zr_3 = 0
            E_96Zr_3 =  0

        if os.path.isfile(u_96Zr_empire4): # Y-91
            #print("Ir 191 file: ",f_191Ir)
            print("u_96Zr_Y-91 exists")
            CS_96Zr_4 = np.genfromtxt(u_96Zr_empire4, delimiter='	', usecols=[1])
            E_96Zr_4 = np.genfromtxt(u_96Zr_empire4, delimiter='	', usecols=[0])
        else:
            print("96Zr_Y-91 file does not exist")
            CS_96Zr_4 = 0
            E_96Zr_4 =  0

        if os.path.isfile(u_96Zr_empire5): # Y-91M
            #print("Ir 191 file: ",f_191Ir)
            print("u_96Zr_Y-91M exists")
            CS_96Zr_5 = np.genfromtxt(u_96Zr_empire5, delimiter='	', usecols=[1])
            E_96Zr_5 = np.genfromtxt(u_96Zr_empire5, delimiter='	', usecols=[0])
        else:
            print("96Zr_Y-91M file does not exist")
            CS_96Zr_5 = 0
            E_96Zr_5 =  0

        if os.path.isfile(u_96Zr_empire6): # Y-92
            #print("Ir 191 file: ",f_191Ir)
            print("u_96Zr_Y-92 exists")
            CS_96Zr_6 = np.genfromtxt(u_96Zr_empire6, delimiter='	', usecols=[1])
            E_96Zr_6 = np.genfromtxt(u_96Zr_empire6, delimiter='	', usecols=[0])
        else:
            print("96Zr_Y-92 file does not exist")
            CS_96Zr_6 = 0
            E_96Zr_6 =  0

        if os.path.isfile(u_96Zr_empire7): # Y-95
            #print("Ir 191 file: ",f_191Ir)
            print("u_96Zr_Y-95 exists")
            CS_96Zr_7 = np.genfromtxt(u_96Zr_empire7, delimiter='	', usecols=[1])
            E_96Zr_7 = np.genfromtxt(u_96Zr_empire7, delimiter='	', usecols=[0])
        else:
            print("96Zr_Y-95 file does not exist")
            CS_96Zr_7 = 0
            E_96Zr_7 =  0

        if os.path.isfile(u_96Zr_empire8): # Zr-95
            #print("Ir 191 file: ",f_191Ir)
            print("u_96Zr_Zr-95 exists")
            CS_96Zr_8 = np.genfromtxt(u_96Zr_empire8, delimiter='	', usecols=[1])
            E_96Zr_8 = np.genfromtxt(u_96Zr_empire8, delimiter='	', usecols=[0])
        else:
            print("96Zr_Zr-95 file does not exist")
            CS_96Zr_8 = 0
            E_96Zr_8 =  0

        E_90Zr = (E_90Zr_1 + E_90Zr_2 + E_90Zr_3) * abund_90Zr
        E_91Zr = (E_91Zr_1 + E_91Zr_2 + E_91Zr_3 + E_91Zr_4 + E_91Zr_5) * abund_91Zr
        E_92Zr = (E_92Zr_1 + E_92Zr_2 + E_92Zr_3 + E_92Zr_4 + E_92Zr_5 + E_92Zr_6) * abund_92Zr
        E_94Zr = (E_94Zr_1 + E_94Zr_2 + E_94Zr_3 + E_94Zr_4 + E_94Zr_5 + E_94Zr_6 + E_94Zr_7) * abund_94Zr
        E_96Zr = (E_96Zr_1 + E_96Zr_2 + E_96Zr_3 + E_96Zr_4 + E_96Zr_5 + E_96Zr_6 + E_96Zr_7 + E_96Zr_8) * abund_96Zr

        CS_90Zr = (CS_90Zr_1 + CS_90Zr_2 + CS_90Zr_3) * abund_90Zr
        CS_91Zr = (CS_91Zr_1 + CS_91Zr_2 + CS_91Zr_3 + CS_91Zr_4 + CS_91Zr_5) * abund_91Zr
        CS_92Zr = (CS_92Zr_1 + CS_92Zr_2 + CS_92Zr_3 + CS_92Zr_4 + CS_92Zr_5 + CS_92Zr_6) * abund_92Zr
        CS_94Zr = (CS_94Zr_1 + CS_94Zr_2 + CS_94Zr_3 + CS_94Zr_4 + CS_94Zr_5 + CS_94Zr_6 + CS_94Zr_7) * abund_94Zr
        CS_96Zr = (CS_96Zr_1 + CS_96Zr_2 + CS_96Zr_3 + CS_96Zr_4 + CS_96Zr_5 + CS_96Zr_6 + CS_96Zr_7 + CS_96Zr_8) * abund_96Zr

        E = E_90Zr + E_91Zr + E_92Zr + E_94Zr + E_96Zr
        CS = CS_90Zr + CS_91Zr + CS_92Zr + CS_94Zr + CS_96Zr

        plt.plot(E, CS, 'y--', label='EMPIRE')





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
    TALYS(foil, A, Z)

    EXFOR(foil, reaction)
    #no EXFOR for :
    # natZn(n,x)67Cu
    Tendl(foil, A, Z)
    #CoH(foil, A, Z, filename_CS)
    EMPIRE(foil, A, Z, filename_CS)

# integral for each of the codes


    mydata(filename_CS, foil)

    plt.legend()
    plt.title(reaction)
    plt.xlim(0,40)
    plt.ylim(ymin=0)
    plt.xlabel('Neutron Energy (MeV)')
    plt.ylabel('Cross Section (mb)')
    plt.savefig('Cross_section_curves/'+ reaction + '.png', dpi=300)
    plt.show()


# LINE 22 - ALICE
# LINE 77 - TALYS
# LINE 100 - EXFOR
# LINE 150 - TENDL
# LINE 441 - CoH
# LINE 489 - EMPIRE
# LINE 942 - Cross_section


#Cross_section('foil', 'A', 'Z', 'EXFOR', 'mydata')
#Cross_section('foil', 'A', 'Z', 'reaction', 'filename_CS')

### ZINK ###

#Cross_section('Zn', '62', '30', 'Zn(n,x)62Zn', '62ZN')
#Cross_section('Zn', '63', '30', 'Zn(n,x)63Zn', '63ZN')
#Cross_section('Zn', '64', '29', 'nat-Zn(n,x)64Cu', '64CU')
#Cross_section('Zn', '65', '28', 'Zn(n,x)65Ni', '65NI')
#Cross_section('Zn', '65', '30', 'Zn(n,x)65Zn', '65ZN')
#Cross_section('Zn', '66', '29', 'Zn(n,x)66Cu', '66CU')
Cross_section('Zn', '67', '29', '67Zn(n,x)67Cu', '67CU') # really 68Zn not naturall in exfor
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

#Cross_section('In', '113', '49', 'In(n,x)111In', '111IN')
#Cross_section('In', '115', '49', '113In(n,2n)112In', '112IN')
#Cross_section('In', '113', '49', 'In(n,x)112mIn', '112INm')
#Cross_section('In', '115', '49', 'In(n,x)112mIn', '112INm')


###  Aluminum ###

#Cross_section('Al', '27', '13', '27Al(n,x)24Na', '24Na')

### Yttrium ###

#Cross_section('Y', '89', '39', '88Y(n,2n)87Y', '87Y')
#Cross_section('Y', '89', '39', '', '87Ym')
#Cross_section('Y', '89', '39', '89Y(n,g)90mY', '90Ym')
#Cross_section('Y', '89', '39', '', '87SRm')
