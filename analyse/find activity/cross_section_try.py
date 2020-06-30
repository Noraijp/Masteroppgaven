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
    print('------------------ ALICE ----------------')
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
        use_cols = 5

        print('ALICE A : ', A)
        print('ALICE Z : ', Z)
        if foil == 'Zn':
            num_header_line = 59
            if A == '69' and Z == '30':
                use_cols = 7
        elif foil == 'Zr':
            num_header_line = 59
            if A == '90' and Z == '39':
                use_cols = 7
            elif A == '91' and Z == '39':
                use_cols = 7
        elif foil == 'Al':
            num_header_line = 51
        elif foil == 'Y':
            num_header_line = 51
        elif foil == 'In':
            num_header_line = 56

        print('use_cols : ', use_cols)

        E  = np.genfromtxt(filename, delimiter=' ', usecols=[0], skip_header=num_header_line, skip_footer=(len(content_full)-len(content)))
        #Z_  = np.genfromtxt(filename, delimiter=' ', usecols=[1], skip_header=56, skip_footer=(len(content_full)-len(content)))
        #A_  = np.genfromtxt(filename, delimiter=' ', usecols=[2], skip_header=56, skip_footer=(len(content_full)-len(content)))
        #print(E)
        CS = np.genfromtxt(filename, delimiter=' ', usecols=[use_cols], skip_header=num_header_line, skip_footer=(len(content_full)-len(content)))
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
        print("E ALICE: ",E_new)
        print("CS ALICE: ", CS_new)
        #plt.plot(E_new, xaxis, label='ALICE')
        tck = interpolate.splrep(E_new, CS_new, s=0)
        E = np.linspace(0, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)



        plt.plot(E, CS, color='m', label='ALICE-2017')
        #plt.show()
        return E, CS, tck




path_ = '/Users/Nora/Documents/UiO/Masteroppgaven/analyse/'


def TALYS(foil, A, Z, file_ending='.tot'):
    # Z = 0XX, A=0XX

    if foil == 'Zn':
        if Z == '030' and A == '069':
            filename = path + '../Talys/' +foil+ '/rp' + Z  + A + '.L01'
        else:
            filename = path + '../Talys/' +foil+ '/rp' + Z  + A + file_ending
        #print('Filename TALYS: ', filename)

    elif foil == 'Al':
    #Al
        filename = path + '../Talys/' +foil+ '/rp'+ Z + A + file_ending + '.027' #'.tot'

    elif  foil == 'In':
    # In
        filename = path + '../Talys/' +foil+ '/rp'+ Z + A + file_ending  #'.tot

    elif foil == 'Y':
    # Y
        filename = path + '../Talys/' +foil+ '/rp'+ Z + A + file_ending + '.089' #'.tot

    #Zr
    elif foil == 'Zr':
        if Z == '039' and A == '090':

        #finn abundance til hvert isotop fra nndc
            abund_90Zr = 0.5145 ; abund_91Zr = 0.1122 ; abund_92Zr = 0.1715 ; abund_94Zr = 0.1738 ; abund_96Zr = 0.280 ;

            print('TALYS Y90m : ')
            f_90Zr = path + '../Tendl/' + foil + '/rp' + Z + A + '.L02' + '.090'
            f_91Zr = path + '../Talys/' + foil + '/rp' + Z  + A + '.L02' + '.091'
            f_92Zr = path + '../Talys/' + foil + '/rp' + Z  + A + '.L02' + '.092'
            f_94Zr = path + '../Talys/' + foil + '/rp' + Z  + A + '.L02' + '.094'
            f_96Zr = path + '../Talys/' + foil + '/rp' + Z  + A + '.L02' + '.096'

            #print(f_90Zr)
            #print(f_61Ni)
            #print(f_62Ni)
            #print(f_64Ni)
            #E = []; CS = []


            if os.path.isfile(f_90Zr):
                #print("Ir 191 file: ",f_191Ir)
                print("f_90Zr exists")
                CS_90Zr = np.genfromtxt(f_90Zr, delimiter=' ', usecols=[1],skip_header=5) * abund_90Zr
                E_Zr = np.genfromtxt(f_90Zr, delimiter=' ', usecols=[0],skip_header=5)
                print(len(E_90Zr))
            else:
                print("90Zr file does not exist")
                CS_90Zr = 0
                E_90Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

            if os.path.isfile(f_91Zr):
                #print("Ir 191 file: ",f_191Ir)
                print("f_91Zr exists")
                CS_91Zr = np.genfromtxt(f_91Zr, delimiter=' ', usecols=[1],skip_header=5) * abund_91Zr
                E_Zr = np.genfromtxt(f_91Zr, delimiter=' ', usecols=[0],skip_header=5)
            else:
                print("91Zr file does not exist")
                CS_91Zr = 0
                E_91Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

            if os.path.isfile(f_92Zr):
                #print("Ir 191 file: ",f_191Ir)
                print("f_92Zr exists")
                CS_92Zr = np.genfromtxt(f_92Zr, delimiter=' ', usecols=[1],skip_header=5) * abund_92Zr
                E_Zr = np.genfromtxt(f_92Zr, delimiter=' ', usecols=[0],skip_header=5)
            else:
                print("92Zr file does not exist")
                CS_92Zr = 0
                E_92Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

            if os.path.isfile(f_94Zr):
                #print("Ir 191 file: ",f_191Ir)
                print("f_94Zr exists")
                CS_94Zr = np.genfromtxt(f_94Zr, delimiter=' ', usecols=[1],skip_header=5) * abund_94Zr
                E_Zr = np.genfromtxt(f_94Zr, delimiter=' ', usecols=[0],skip_header=5)
            else:
                print("94Zr file does not exist")
                CS_94Zr = 0
                E_94Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

            if os.path.isfile(f_96Zr):
                #print("Ir 191 file: ",f_191Ir)
                print("f_96Zr exists")
                CS_96Zr = np.genfromtxt(f_96Zr, delimiter=' ', usecols=[1],skip_header=5) * abund_96Zr
                E_Zr = np.genfromtxt(f_96Zr, delimiter=' ', usecols=[0],skip_header=5)
            else:
                print("96Zr file does not exist")
                CS_96Zr = 0
                E_96Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

            #print(E_58Ni)
            #CS_193Ir = np.genfromtxt(f_193Ir, delimiter=' ', usecols=[1],skip_header=5)
            #E_193Ir = np.genfromtxt(f_193Ir, delimiter=' ', usecols=[0],skip_header=5)

            Ee = E_Zr
            #CS = CS_90Zr*abund_90Zr + CS_91Zr*abund_91Zr + CS_92Zr*abund_92Zr + CS_94Zr*abund_94Zr + CS_96Zr*abund_96Zr
            CSs = CS_90Zr + CS_91Zr + CS_92Zr + CS_94Zr + CS_96Zr

            #print('E_zr :', len(E))
            #print('CS_zr :', len(CS))
            #CS = CS_191Ir*abund_191Ir + CS_193Ir*abund_193Ir

            #plt.plot(E, xaxis, 'y--', label='TENDL')

            # tck = interpolate.splrep(E_, CS_, s=0)
            # E = np.linspace(0, 40, 1000)
            # CS = interpolate.splev(E, tck, der=0)
            # plt.plot(E, CS, color='r', label='TENDL-2019')
            # #plt.plot(E_61Ni,CS_61Ni, label='61Ni', linewidth=0.5)
            #
            # return E, CS, tck

        elif Z == '039' and A == '091':
            filename = path + '../Talys/' +foil+ '/rp' + Z  + A + '.L01'

            Ee  = np.genfromtxt(filename, delimiter=' ', usecols=[0],skip_header=5)
            CSs = np.genfromtxt(filename, delimiter=' ', usecols=[1],skip_header=5)

            print('TALYS filename :', filename)

        else:
            filename = path + '../Talys/' +foil+ '/rp'+ Z + A + file_ending #'.tot

            Ee  = np.genfromtxt(filename, delimiter=' ', usecols=[0],skip_header=5)
            CSs = np.genfromtxt(filename, delimiter=' ', usecols=[1],skip_header=5)

            print('TALYS filename :', filename)

    #filename = self.path + '/../Talys/' +foil+ '/rp'+Z+A+'.L02'
    #Ee  = np.genfromtxt(filename, delimiter=' ', usecols=[0],skip_header=5)
    #CSs = np.genfromtxt(filename, delimiter=' ', usecols=[1],skip_header=5)


    #print(CS)
    print('TALYS CS: ', CSs)
    print(Ee)

    tck = interpolate.splrep(Ee, CSs, s=0)
    E = np.linspace(0, 40, 1000)
    CS = interpolate.splev(E, tck, der=0)
    plt.plot(E, CS, color='g', label='TALYS-1.9')

    #plt.show()
    return E, CS, tck



def EXFOR(foil, reaction):
    filename = path + '/../EXFOR/' + foil +'/' + reaction + '.txt'
    #print(filename)
    print(reaction)

    if os.path.isfile(filename):

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
            plt.errorbar(E, CS, marker='.', markersize=4, linewidth=0.001, xerr=dE, yerr=dCS, elinewidth=0.5, capthick=0.5, capsize=3.0, label=author[0] )
            #plt.legend()
            #plt.show()

            return E, dE, CS, dCS, author
        # else:
        #     print("exfor file does not exist for {}".format(reaction))
        #     return 0, 0, 0, 0, '0'


def Tendl(foil, A, Z, file_ending='.tot'):
    num_files = 0
    print('TENDL!!!!!!!!!!!!!!!')
    print("foil: ",foil )
    print("Z: ", Z )
    print("A: ", A  )

    if foil == 'Zn':
        #A = ['191', '193'] # stable iridium isotopes
        abund_64Zn = 0.4917 ; abund_66Zn = 0.2773 ; abund_67Zn = 0.404 ; abund_68Zn = 0.1845 ; abund_70Zn = 0.061;
        #file_ending =
        #endre f_191I til de stabile isotopene i zink osv

        if Z == '030' and A == '069':
            print('TENDL zn69m')
            f_64Zn = path + '/../Tendl/' + foil + '/rp030064_' + Z + A + '.L01' + '.txt'
            f_66Zn = path + '/../Tendl/' + foil + '/rp030066_' + Z + A + '.L01' + '.txt'
            f_67Zn = path + '/../Tendl/' + foil + '/rp030067_' + Z + A + '.L01' + '.txt'
            f_68Zn = path + '/../Tendl/' + foil + '/rp030068_' + Z + A + '.L01' + '.txt'
            f_70Zn = path + '/../Tendl/' + foil + '/rp030070_' + Z + A + '.L01' + '.txt'
        else:
            f_64Zn = path + '/../Tendl/' + foil + '/rp030064_' + Z + A + file_ending + '.txt'
            f_66Zn = path + '/../Tendl/' + foil + '/rp030066_' + Z + A + file_ending + '.txt'
            f_67Zn = path + '/../Tendl/' + foil + '/rp030067_' + Z + A + file_ending + '.txt'
            f_68Zn = path + '/../Tendl/' + foil + '/rp030068_' + Z + A + file_ending + '.txt'
            f_70Zn = path + '/../Tendl/' + foil + '/rp030070_' + Z + A + file_ending + '.txt'
        print('Filename TENDL: ', f_70Zn)


        #print("Ir 193 file: ",f_193Ir)
        #print("Ir 191 file: ",f_191Ir)
        if os.path.isfile(f_64Zn):
            #print("Ir 191 file: ",f_191Ir)
            print("f_64Zn exists")
            CS_64Zn = np.genfromtxt(f_64Zn, delimiter=' ', usecols=[1],skip_header=5)
            E_Zn = np.genfromtxt(f_64Zn, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("f_64Zn file does not exist")
            CS_64Zn = 0
            E_64Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        if os.path.isfile(f_66Zn):
            #print("Ir 191 file: ",f_191Ir)
            #print("f_191Ir exists")
            CS_66Zn = np.genfromtxt(f_66Zn, delimiter=' ', usecols=[1],skip_header=5)
            E_Zn = np.genfromtxt(f_66Zn, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("Zn 66 file does not exist")
            CS_66Zn = 0
            E_66Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        if os.path.isfile(f_67Zn):
            #print("Ir 191 file: ",f_191Ir)
            print(" Zn 67 exists")
            CS_67Zn = np.genfromtxt(f_67Zn, delimiter=' ', usecols=[1],skip_header=5)
            E_Zn = np.genfromtxt(f_67Zn, delimiter=' ', usecols=[0],skip_header=5)
            print(CS_67Zn)
        else:
            print("f_67Zn file does not exist")
            CS_67Zn = 0
            E_67Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        #print('CS_67Zn file :', CS_67Zn)


        if os.path.isfile(f_68Zn):
            #print("Ir 191 file: ",f_191Ir)
            print("f_68Zn exists")
            CS_68Zn = np.genfromtxt(f_68Zn, delimiter=' ', usecols=[1],skip_header=5)
            E_Zn = np.genfromtxt(f_68Zn, delimiter=' ', usecols=[0],skip_header=5)
            #print(CS_68Zn)
        else:
            print("f_68Zn file does not exist")
            CS_68Zn = 0
            E_68Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        #print('CS_68Zn file :', CS_68Zn)

        if os.path.isfile(f_70Zn):
            #print("Ir 191 file: ",f_191Ir)
            print("f_70Zn exists")
            CS_70Zn = np.genfromtxt(f_70Zn, delimiter=' ', usecols=[1],skip_header=5)
            E_Zn = np.genfromtxt(f_70Zn, delimiter=' ', usecols=[0],skip_header=5)
            print(f_70Zn)
        else:
            print("f_70Zn file does not exist")
            CS_70Zn = 0
            E_70Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        #print('CS_70Zn file :', CS_70Zn)



        print('Number of files : ' ,num_files)

        En_Zn = E_Zn
        CS_Zn = CS_64Zn*abund_64Zn + CS_66Zn*abund_66Zn + CS_67Zn*abund_67Zn + CS_68Zn*abund_68Zn + CS_70Zn*abund_70Zn

        #print('Energy TENDL :', E_)
        #print('CS TENDL :', CS_)
        print('...........................')
        tck = interpolate.splrep(En_Zn, CS_Zn, s=0)
        E = np.linspace(1, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)
        plt.plot(E, CS, color='r', label='TENDL-2019')

        #plt.plot(E_191Ir,CS_191Ir, label='191Ir')
        #plt.plot(E_193Ir,CS_193Ir, label='193Ir')
        #plt.plot(E, CS, label='tot')
        #plt.legend()
        #plt.show()
        return E, CS, tck



    elif foil == 'Zr':
        #finn abundance til hvert isotop fra nndc
        abund_90Zr = 0.5145 ; abund_91Zr = 0.1122 ; abund_92Zr = 0.1715 ; abund_94Zr = 0.1738 ; abund_96Zr = 0.280 ;

        if Z == '039' and A == '090':
            print('TENDL Y90m : ')
            f_90Zr = path + '/../Tendl/' + foil + '/rp040090_' + Z + A + '.L02' + '.txt'
            f_91Zr = path + '/../Tendl/' + foil + '/rp040091_' + Z  + A + '.L02' + '.txt'
            f_92Zr = path + '/../Tendl/' + foil + '/rp040092_' + Z  + A + '.L02' + '.txt'
            f_94Zr = path + '/../Tendl/' + foil + '/rp040094_' + Z  + A + '.L02' + '.txt'
            f_96Zr = path + '/../Tendl/' + foil + '/rp040096_' + Z  + A + '.L02' + '.txt'

        if Z == '039' and A == '091':
            print('TENDL Y91m : ')
            f_90Zr = path + '/../Tendl/' + foil + '/rp040090_' + Z + A + '.L01' + '.txt'
            f_91Zr = path + '/../Tendl/' + foil + '/rp040091_' + Z  + A + '.L01' + '.txt'
            f_92Zr = path + '/../Tendl/' + foil + '/rp040092_' + Z  + A + '.L01' + '.txt'
            f_94Zr = path + '/../Tendl/' + foil + '/rp040094_' + Z  + A + '.L01' + '.txt'
            f_96Zr = path + '/../Tendl/' + foil + '/rp040096_' + Z  + A + '.L01' + '.txt'

        else:
        #f_90Zr = path + '/../Tendl/' + foil + '/rp040090_' + Z  + A + file_ending + '.txt'
            f_90Zr = path + '/../Tendl/' + foil + '/rp040090_' + Z + A + file_ending + '.txt'
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
            E_Zr = np.genfromtxt(f_90Zr, delimiter=' ', usecols=[0],skip_header=5)
            print('CS_90Zr : ', len(CS_90Zr))

        else:
            print("90Zr file does not exist")
            CS_90Zr = 0
            E_90Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        if os.path.isfile(f_91Zr):
            #print("Ir 191 file: ",f_191Ir)
            print("f_91Zr exists")
            CS_91Zr = np.genfromtxt(f_91Zr, delimiter=' ', usecols=[1],skip_header=5) * abund_91Zr
            E_Zr = np.genfromtxt(f_91Zr, delimiter=' ', usecols=[0],skip_header=5)
            print('CS_91Zr : ', len(CS_91Zr))
        else:
            print("91Zr file does not exist")
            CS_91Zr = 0
            E_91Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        if os.path.isfile(f_92Zr):
            #print("Ir 191 file: ",f_191Ir)
            print("f_92Zr exists")
            CS_92Zr = np.genfromtxt(f_92Zr, delimiter=' ', usecols=[1],skip_header=5) * abund_92Zr
            E_Zr = np.genfromtxt(f_92Zr, delimiter=' ', usecols=[0],skip_header=5)
            print('CS_92Zr : ', len(CS_92Zr))
        else:
            print("92Zr file does not exist")
            CS_92Zr = 0
            E_92Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        if os.path.isfile(f_94Zr):
            #print("Ir 191 file: ",f_191Ir)
            print("f_94Zr exists")
            CS_94Zr = np.genfromtxt(f_94Zr, delimiter=' ', usecols=[1],skip_header=5) * abund_94Zr
            E_Zr = np.genfromtxt(f_94Zr, delimiter=' ', usecols=[0],skip_header=5)
            print('CS_94Zr : ', len(CS_94Zr))
        else:
            print("94Zr file does not exist")
            CS_94Zr = 0
            E_94Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        if os.path.isfile(f_96Zr):
            #print("Ir 191 file: ",f_191Ir)
            print("f_96Zr exists")
            CS_96Zr = np.genfromtxt(f_96Zr, delimiter=' ', usecols=[1],skip_header=5) * abund_96Zr
            E_Zr = np.genfromtxt(f_96Zr, delimiter=' ', usecols=[0],skip_header=5)
            print('CS_96Zr : ', len(CS_96Zr))
        else:
            print("96Zr file does not exist")
            CS_96Zr = 0
            E_96Zr =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        #print(E_58Ni)
        #CS_193Ir = np.genfromtxt(f_193Ir, delimiter=' ', usecols=[1],skip_header=5)
        #E_193Ir = np.genfromtxt(f_193Ir, delimiter=' ', usecols=[0],skip_header=5)

        E_ = E_Zr
        #CS = CS_90Zr*abund_90Zr + CS_91Zr*abund_91Zr + CS_92Zr*abund_92Zr + CS_94Zr*abund_94Zr + CS_96Zr*abund_96Zr
        CS_ = CS_90Zr + CS_91Zr + CS_92Zr + CS_94Zr + CS_96Zr

        #print('E_zr :', len(E))
        #print('CS_zr :', len(CS))
        #CS = CS_191Ir*abund_191Ir + CS_193Ir*abund_193Ir

        #plt.plot(E, xaxis, 'y--', label='TENDL')

        tck = interpolate.splrep(E_, CS_, s=0)
        E = np.linspace(0, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)
        plt.plot(E, CS, color='r', label='TENDL-2019')
        #plt.plot(E_61Ni,CS_61Ni, label='61Ni', linewidth=0.5)

        return E, CS, tck



    elif foil == 'Y':
        #finn abundance til hvert isotop fra nndc
        abund_89Y = 1.00;


        f_89Y = path + '/../Tendl/' + foil + '/rp039089_' + Z + A + file_ending + '.txt'

        print("Y 89 file: ",f_89Y)
        if os.path.isfile(f_89Y):
            #print("Ir 191 file: ",f_191Ir)
            print("f_89Y exists")
            CS_89Y = np.genfromtxt(f_89Y, delimiter=' ', usecols=[1],skip_header=5)
            E_89Y = np.genfromtxt(f_89Y, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("89Y file does not exist")
            CS_89Y = 0
            E_89Y =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        E_ = E_89Y
        CS_ = CS_89Y*abund_89Y

        tck = interpolate.splrep(E_, CS_, s=0)
        E = np.linspace(0, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)
        plt.plot(E, CS, color='r', label='TENDL-2019')

        return E, CS, tck


    elif foil == 'In':
        #finn abundance til hvert isotop fra nndc
        abund_113In = 0.429;  abund_115In = 0.9571;


        f_113In = path + '/../Tendl/' + foil + '/rp049113_' + Z + A + file_ending + '.txt'
        f_115In = path + '/../Tendl/' + foil + '/rp049115_' + Z + A + file_ending + '.txt'

        print("In 113 file: ",f_113In)
        if os.path.isfile(f_113In):
            #print("Ir 191 file: ",f_191Ir)
            print("f_113In exists")
            CS_113In = np.genfromtxt(f_113In, delimiter=' ', usecols=[1],skip_header=5)
            E_In = np.genfromtxt(f_113In, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("f_113In file does not exist")
            CS_113In = 0
            E_113In =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

        print("In 115 file: ",f_115In)
        if os.path.isfile(f_115In):
            #print("Ir 191 file: ",f_191Ir)
            print("f_115In exists")
            CS_115In = np.genfromtxt(f_115In, delimiter=' ', usecols=[1],skip_header=5)
            E_In = np.genfromtxt(f_115In, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("f_115In file does not exist")
            CS_115In = 0
            E_115In =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        CS_ = CS_113In*abund_113In + CS_115In*abund_115In

        tck = interpolate.splrep(E_In, CS_, s=0)
        E = np.linspace(0, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)
        plt.plot(E, CS, color='r', label='TENDL-2019')


        return E, CS, tck


    elif foil == 'Al':
        #finn abundance til hvert isotop fra nndc
        abund_27Al = 1.00;


        f_27Al = path + '/../Tendl/' + foil + '/rp013027_' + Z + A + file_ending + '.txt'

        print('f_27AL, TENDL :', f_27Al)
        if os.path.isfile(f_27Al):
            #print("Ir 191 file: ",f_191Ir)
            print("f_27Al exists")
            CS_27Al = np.genfromtxt(f_27Al, delimiter=' ', usecols=[1],skip_header=5)
            E_27Al = np.genfromtxt(f_27Al, delimiter=' ', usecols=[0],skip_header=5)
        else:
            print("27Al file does not exist, TENDL")
            CS_27Al = 0
            E_27Al =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


        E_ = E_27Al
        CS_ = CS_27Al*abund_27Al

        #plt.plot(E, xaxis, 'y--', label='TENDL')

        tck = interpolate.splrep(E_, CS_, s=0)
        E = np.linspace(0, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)
        plt.plot(E, CS, color='r', label='TENDL-2019')

        return E, CS, tck





def CoH(foil, A, Z, filename_CS_me):
    print('----------CoH----------')
    #print("foil: ",foil )
    #print("Z: ", Z )
    #print("A: ", A  )

    if foil == 'Zn':
        #A = ['191', '193'] # stable iridium isotopes
        abund_64Zn = 0.4917 ; abund_66Zn = 0.2773 ; abund_67Zn = 0.404 ; abund_68Zn = 0.1845 ; abund_70Zn = 0.061;
        #file_ending =
        #endre f_191I til de stabile isotopene i zink osv
        v_64Zn_Coh = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-0' + filename_CS_me + '_coh' + '.txt'
        v_66Zn_Coh = path + '/../EMPIRECOH2/' + foil + '/66Zn/' +  Z + '-0' + filename_CS_me + '_coh' + '.txt'
        u_67Zn_Coh = path + '/../EMPIRECOH2/' + foil + '/67Zn/' +  Z + '-0' + filename_CS_me + '_coh' + '.txt'
        u_68Zn_Coh = path + '/../EMPIRECOH2/' + foil + '/68Zn/' +  Z + '-0' + filename_CS_me + '_coh' + '.txt'


        #v_64Zn_empire = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-' + 'Cu-64' + '_empire' + '.txt'
        v_67Zn = path + '/../EMPIRECOH2/' + foil + '/67Zn/' +  Z + '-0' + filename_CS_me + '_coh' + '.txt'        #f_67Zn = path + '/../EMPIRECOH2/' + foil + '/rp030067_' + Z +A + file_ending + '.txt'
        #f_68Zn = path + '/../EMPIRECOH2/' + foil + '/rp030068_' + Z +A + file_ending + '.txt'
        #f_70Zn = path + '/../EMPIRECOH2/' + foil + '/rp030070_' + Z +A + file_ending + '.txt'

################### 64Zn #########################

        #print("Zn 64 file: ",v_64Zn_Coh)
        #print("Ir 191 file: ",f_191Ir)
        if os.path.isfile(v_64Zn_Coh):
            #print("Ir 191 file: ",f_191Ir)
            print("v_64Zn exists")
            CS_64Zn_Coh = np.genfromtxt(v_64Zn_Coh, delimiter='\t', usecols=[1])
            E_Zn_Coh = np.genfromtxt(v_64Zn_Coh, delimiter='\t', usecols=[0])
        else:
            print("64 Zn file does not exist")
            CS_64Zn_Coh = 0
            E_64Zn_Coh =  0 #np.genfromtxt(f_191Ir, '\t', usecols=[0],skip_header=5)
        #print('!!!!!!!!!!!!!!!!!!!!!!!!!')
        #print('E_64Zn_Coh :', E_Zn_Coh)
        #print('CS_64Zn_Coh :', CS_64Zn_Coh)
        #plt.plot(E_64Zn_Coh, CS_64Zn_Coh, color='c', label='CoH')

################### 66Zn #########################

        if os.path.isfile(v_66Zn_Coh): # Cu-64
            #print("Ir 191 file: ",f_191Ir)
            print("v_66Zn_Cu-64 exists")
            CS_66Zn_Coh = np.genfromtxt(v_66Zn_Coh, delimiter='\t', usecols=[1])
            E_Zn_Coh = np.genfromtxt(v_66Zn_Coh, delimiter='\t', usecols=[0])
        else:
            print("66Zn_Cu-64 file does not exist")
            CS_66Zn_Coh = 0
            E_66Zn_Coh =  0 #np.genfromtxt(f_191Ir, delimiter='\t', usecols=[0],skip_header=5)

################### 67Zn #########################

        if os.path.isfile(u_67Zn_Coh): # Cu-64
            #print("Ir 191 file: ",f_191Ir)
            print("u_67Zn_Cu-64 exists")
            CS_67Zn_Coh = np.genfromtxt(u_67Zn_Coh, delimiter='\t', usecols=[1])
            E_Zn_Coh = np.genfromtxt(u_67Zn_Coh, delimiter='\t', usecols=[0])
        else:
            print("67Zn_Cu-64 file does not exist")
            CS_67Zn_Coh = 0
            E_67Zn_Coh =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


################### 68Zn #########################

        if os.path.isfile(u_68Zn_Coh): # Cu-64
            #print("Ir 191 file: ",f_191Ir)
            print("u_68Zn_Cu-64 exists")
            CS_68Zn_Coh = np.genfromtxt(u_68Zn_Coh, delimiter='\t', usecols=[1])
            E_Zn_Coh = np.genfromtxt(u_68Zn_Coh, delimiter='\t', usecols=[0])
        else:
            print("68Zn_Cu-64 file does not exist")
            CS_68Zn_Coh = 0
            E_68Zn_Coh =  0


        #plt.plot(E_67Zn_Coh, CS_67Zn_Coh, color='c', label='CoH')
        CS_Zn = CS_64Zn_Coh*abund_64Zn + CS_66Zn_Coh*abund_66Zn + CS_67Zn_Coh*abund_67Zn + CS_68Zn_Coh*abund_68Zn
        E_Zn = E_Zn_Coh
        #print('CS_64Cu Coh :', CS_64Cu)
        #print('E_64Cu Coh :', E_64Cu)

        # to make the plot smoother
        tck = interpolate.splrep(E_Zn, CS_Zn, s=0) #interpalation - give it raw data of energy array and cross section array that we have calculated
        E = np.linspace(0, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)
        plt.plot(E, CS, color='c', label='CoH-3.5.3')

        return E, CS, tck

#-------------------------- Al -------------------------------

    if foil == 'Al':
        abund_27Al = 1
        v_27Al = path + '/../EMPIRECOH2/' + foil + '/27Al/' +  Z + '-0' + filename_CS_me + '_coh' + '.txt'

        if os.path.isfile(v_27Al):
            #print("Ir 191 file: ",f_191Ir)
            print("v_27Al exists")
            CS_27Al = np.genfromtxt(v_27Al, delimiter='\t', usecols=[1])
            E_27Al = np.genfromtxt(v_27Al, delimiter='\t', usecols=[0])
        else:
            print("27Al file does not exist")
            CS_27Al = 0
            E_27Al =  0 #np.genfromtxt(f_191Ir, delimiter='\t', usecols=[0],skip_header=5)


        CS_Al = CS_27Al*abund_27Al

        # to make the plot smoother
        tck = interpolate.splrep(E_27Al, CS_Al, s=0) #interpalation - give it raw data of energy array and cross section array that we have calculated
        E = np.linspace(0, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)
        plt.plot(E, CS, color='c', label='CoH-3.5.3')

        return E, CS, tck

#----------------------------- In ------------------------

    if foil == 'In':
        abund_113In = 0.429 ; abund_115In = 0.9571

        v_113In_Coh = path + '/../EMPIRECOH2/' + foil + '/113In/' +  Z + '-' + filename_CS_me + '_coh' + '.txt'
        v_115In_Coh = path + '/../EMPIRECOH2/' + foil + '/115In/' +  Z + '-' + filename_CS_me + '_coh' + '.txt'

        print("In 113 file: ",v_113In_Coh)
        if os.path.isfile(v_113In_Coh):
            #print("Ir 191 file: ",f_191Ir)
            #print("v_113In_Coh exists")
            CS_113In_Coh = np.genfromtxt(v_113In_Coh, delimiter='\t', usecols=[1])
            E_In_Coh = np.genfromtxt(v_113In_Coh, delimiter='\t', usecols=[0])
        else:
            #print("27Al file does not exist")
            CS_113In_Coh = 0
            E_113In_Coh =  0 #np.genfromtxt(f_191Ir, delimiter='\t', usecols=[0],skip_header=5)

        print("In 115 file: ",v_115In_Coh)
        if os.path.isfile(v_115In_Coh):
            #print("Ir 191 file: ",f_191Ir)
            #print("v_115In_Coh exists")
            CS_115In_Coh = np.genfromtxt(v_115In_Coh, delimiter='\t', usecols=[1])
            E_In_Coh = np.genfromtxt(v_115In_Coh, delimiter='\t', usecols=[0])
        else:
            print("27Al file does not exist")
            CS_115In_Coh = 0
            E_115In_Coh =  0 #np.genfromtxt(f_191Ir, delimiter='\t', usecols=[0],skip_header=5)


        CS_In = CS_113In_Coh*abund_113In + CS_115In_Coh*abund_115In

        # to make the plot smoother
        tck = interpolate.splrep(E_In_Coh, CS_In, s=0) #interpalation - give it raw data of energy array and cross section array that we have calculated
        E = np.linspace(0, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)
        plt.plot(E, CS, color='c', label='CoH-3.5.3')

        return E, CS, tck

#----------------------- Y -----------------------

    if foil == 'Y':
        #A = ['191', '193'] # stable iridium isotopes
        abund_89Y = 1.0
        #file_ending =
        #endre f_191I til de stabile isotopene i zink osv
        #v_64Zn = path + '/../EMPIRECOH2/' + foil + '/folder/' +  Z + '-0' + filename_CS + '_coh' + '.txt'
        # 64Zn
        v_89Y_Coh = path + '/../EMPIRECOH2/' + foil + '/89Y/' +  Z + '-0' + filename_CS_me + '_coh' + '.txt'

        if os.path.isfile(v_89Y_Coh):
            #print("Ir 191 file: ",f_191Ir)
            #print("v_115In_Coh exists")
            CS_89Y_Coh = np.genfromtxt(v_89Y_Coh, delimiter='\t', usecols=[1])
            E_Y_Coh = np.genfromtxt(v_89Y_Coh, delimiter='\t', usecols=[0])
        else:
            print("27Al file does not exist")
            CS_89Y_Coh = 0
            E_89Y_Coh =  0 #np.genfromtxt(f_191Ir, delimiter='\t', usecols=[0],skip_header=5)


        E_Y = E_Y_Coh
        CS_Y = CS_89Y_Coh*abund_89Y

        # to make the plot smoother
        tck = interpolate.splrep(E_Y, CS_Y, s=0) #interpalation - give it raw data of energy array and cross section array that we have calculated
        E = np.linspace(0, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)
        plt.plot(E, CS, color='c', label='CoH-3.5.3')
        return E, CS, tck



#----------------------- Zr  -----------------------


    if foil == 'Zr':
        #A = ['191', '193'] # stable iridium isotopes
        abund_90Zr = 0.5145 ; abund_91Zr = 0.1122 ; abund_92Zr = 0.1715 ; abund_94Zr = 0.1738 ; abund_96Zr = 0.280

        v_90Zr_Coh = path + '/../EMPIRECOH2/' + foil + '/90Zr/' +  Z + '-0' + filename_CS_me + '_coh' + '.txt'
        v_91Zr_Coh = path + '/../EMPIRECOH2/' + foil + '/91Zr/' +  Z + '-0' + filename_CS_me + '_coh' + '.txt'
        v_92Zr_Coh = path + '/../EMPIRECOH2/' + foil + '/92Zr/' +  Z + '-0' + filename_CS_me + '_coh' + '.txt'
        v_94Zr_Coh = path + '/../EMPIRECOH2/' + foil + '/94Zr/' +  Z + '-0' + filename_CS_me + '_coh' + '.txt'
        v_96Zr_Coh = path + '/../EMPIRECOH2/' + foil + '/96Zr/' +  Z + '-0' + filename_CS_me + '_coh' + '.txt'

######################## 90Zr #######################
        print("Zr 90 file: ",v_90Zr_Coh)
        if os.path.isfile(v_90Zr_Coh):
            #print("Zr 90 file: ",v_90Zr_Coh)
            print("v_90Zr_Coh exists")
            CS_90Zr_Coh = np.genfromtxt(v_90Zr_Coh, delimiter='\t', usecols=[1])
            E_Zr_Coh = np.genfromtxt(v_90Zr_Coh, delimiter='\t', usecols=[0])
        else:
            print("90Zr file does not exist")
            CS_90Zr_Coh = 0
            E_90Zr_Coh =  0 #np.genfromtxt(f_191Ir, delimiter='\t', usecols=[0],skip_header=5)

###################### 91Zr #########################

        print("Zr 90 file: ",v_91Zr_Coh)
        if os.path.isfile(v_91Zr_Coh):
            print("Zr 91 file: ",v_91Zr_Coh)
            print("v_91Zr_Coh exists")
            CS_91Zr_Coh = np.genfromtxt(v_91Zr_Coh, delimiter='\t', usecols=[1])
            E_Zr_Coh = np.genfromtxt(v_91Zr_Coh, delimiter='\t', usecols=[0])
        else:
            print("91Zr file does not exist")
            CS_91Zr_Coh = 0
            E_91Zr_Coh =  0 #np.genfromtxt(f_191Ir, delimiter='\t', usecols=[0],skip_header=5)

###################### 92Zr #########################

        print("Zr 90 file: ",v_92Zr_Coh)
        if os.path.isfile(v_92Zr_Coh):
            print("Zr 92 file: ",v_92Zr_Coh)
            print("v_92Zr_Coh exists")
            CS_92Zr_Coh = np.genfromtxt(v_92Zr_Coh, delimiter='\t', usecols=[1])
            E_Zr_Coh = np.genfromtxt(v_92Zr_Coh, delimiter='\t', usecols=[0])
        else:
            print("92Zr file does not exist")
            CS_92Zr_Coh = 0
            E_92Zr_Coh =  0 #np.genfromtxt(f_191Ir, delimiter='\t', usecols=[0],skip_header=5)

###################### 94Zr #########################

        print("Zr 90 file: ",v_94Zr_Coh)
        if os.path.isfile(v_94Zr_Coh):
            print("Zr 94 file: ",v_94Zr_Coh)
            print("v_94Zr_Coh exists")
            CS_94Zr_Coh = np.genfromtxt(v_94Zr_Coh, delimiter='\t', usecols=[1])
            E_Zr_Coh = np.genfromtxt(v_94Zr_Coh, delimiter='\t', usecols=[0])
        else:
            print("94Zr file does not exist")
            CS_94Zr_Coh = 0
            E_94Zr_Coh =  0 #np.genfromtxt(f_191Ir, delimiter='\t', usecols=[0],skip_header=5)

###################### 96Zr #########################

        print("Zr 90 file: ",v_96Zr_Coh)
        if os.path.isfile(v_96Zr_Coh):
            print("Zr 96 file: ",v_96Zr_Coh)
            print("v_96Zr_Coh exists")
            CS_96Zr_Coh = np.genfromtxt(v_96Zr_Coh, delimiter='\t', usecols=[1])
            E_Zr_Coh = np.genfromtxt(v_96Zr_Coh, delimiter='\t', usecols=[0])
        else:
            print("96Zr file does not exist")
            CS_96Zr_Coh = 0
            E_96Zr_Coh =  0 #np.genfromtxt(f_191Ir, delimiter='\t', usecols=[0],skip_header=5)



        CS_Zr = CS_90Zr_Coh*abund_90Zr + CS_91Zr_Coh*abund_91Zr + CS_92Zr_Coh*abund_92Zr + CS_94Zr_Coh*abund_94Zr + CS_96Zr_Coh*abund_96Zr

        # to make the plot smoother
        tck = interpolate.splrep(E_Zr_Coh, CS_Zr, s=0) #interpalation - give it raw data of energy array and cross section array that we have calculated
        E = np.linspace(0, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)
        plt.plot(E, CS, color='c', label='CoH-3.5.3')

        return E, CS, tck



def EMPIRE(foil, A, Z, filename_CS):
    print('----------EMPIRE----------')

    #print("folder: ",folder )
    #print("Z: ", Z )
    #print("A: ", A  )

    if foil == 'Zn':
        #A = ['191', '193'] # stable iridium isotopes
        abund_64Zn = 0.4917 ; abund_66Zn = 0.2773 ; abund_67Zn = 0.404 ; abund_68Zn = 0.1845 ; abund_70Zn = 0.061;
        #file_ending =
        #endre f_191I til de stabile isotopene i zink osv
        #v_64Zn = path + '/../EMPIRECOH2/' + foil + '/folder/' +  Z + '-0' + filename_CS + '_coh' + '.txt'
        # 64Zn
        u_64Zn_empire = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-' + filename_CS + '_empire' + '.txt'
        # 66 Zn
        u_66Zn_empire = path + '/../EMPIRECOH2/' + foil + '/66Zn/' +  Z + '-' + filename_CS + '_empire' + '.txt'
        # 67Zn
        u_67Zn_empire = path + '/../EMPIRECOH2/' + foil + '/67Zn/' +  Z + '-' + filename_CS + '_empire' + '.txt'
        # 86Zn
        u_68Zn_empire = path + '/../EMPIRECOH2/' + foil + '/68Zn/' +  Z + '-' + filename_CS + '_empire' + '.txt'
        # 70Zn
        u_70Zn_empire = path + '/../EMPIRECOH2/' + foil + '/70Zn/' +  Z + '-' + filename_CS + '_empire' + '.txt'


        #print("Zn 64 file: ",u_64Zn_empire1)
        #print("Ir 191 file: ",f_191Ir)

############################ 64Zn #################################

        if os.path.isfile(u_64Zn_empire):
            #print("Ir 191 file: ",f_191Ir)
            print("v_64Zn exists")
            CS_64Zn = np.genfromtxt(u_64Zn_empire, delimiter='\t', usecols=[1])
            E_Zn = np.genfromtxt(u_64Zn_empire, delimiter='\t', usecols=[0])
        else:
            print("64Zn file does not exist")
            CS_64Zn = 0
            E_64Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        #print('!!!!!!!!!!!!!!!!!!!!!!!!!')
        #print(E_64Zn_1)
        #plt.plot(E_64Zn_1, CS_64Zn_1, color='m', label='EMPIRE')
        #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)

############################ 66Zn #################################


        if os.path.isfile(u_66Zn_empire):
            #print("Ir 191 file: ",f_191Ir)
            print("v_66Zn exists")
            CS_66Zn = np.genfromtxt(u_66Zn_empire, delimiter='\t', usecols=[1])
            E_Zn = np.genfromtxt(u_66Zn_empire, delimiter='\t', usecols=[0])
        else:
            print("66Zn file does not exist")
            CS_66Zn = 0
            E_66Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)


############################ 67Zn #################################

        if os.path.isfile(u_67Zn_empire): # Cu-67
            #print("Ir 191 file: ",f_191Ir)
            print("u_67Zn exists")
            CS_67Zn = np.genfromtxt(u_67Zn_empire, delimiter='\t', usecols=[1])
            E_Zn = np.genfromtxt(u_67Zn_empire, delimiter='\t', usecols=[0])
        else:
            print("67Zn file does not exist")
            CS_67Zn = 0
            E_67Zn =  0 #np.genfromtxt(f_191Ir, delimiter=' ', usecols=[0],skip_header=5)
        #print('!!!!!!!!!!!!!!!!!!!!!!!!!')
        #print(E_68Zn_E)
        #plt.plot(E_68Zn_E, CS_68Zn_E, color='m', label='EMPIRE')



############################ 68Zn #################################

        if os.path.isfile(u_68Zn_empire): # Ni-65
            #print("Ir 191 file: ",f_191Ir)
            print("u_68Zn exists")
            CS_68Zn = np.genfromtxt(u_68Zn_empire, delimiter='\t', usecols=[1])
            E_Zn = np.genfromtxt(u_68Zn_empire, delimiter='\t', usecols=[0])
        else:
            print("68Zn file does not exist")
            CS_68Zn = 0
            E_68Zn =  0


############################ 70Zn #################################

        if os.path.isfile(u_70Zn_empire): # Ni-65
            #print("Ir 191 file: ",f_191Ir)
            print("u_70Zn exists")
            CS_70Zn = np.genfromtxt(u_70Zn_empire, delimiter='\t', usecols=[1])
            E_Zn = np.genfromtxt(u_70Zn_empire, delimiter='\t', usecols=[0])
        else:
            print("70Zn file does not exist")
            CS_70Zn = 0
            E_70Zn =  0



        CS_Zn = CS_64Zn*abund_64Zn + CS_66Zn*abund_66Zn + CS_67Zn*abund_67Zn + CS_68Zn*abund_68Zn + CS_70Zn*abund_70Zn
        #E_Zn = E_Zn

        #CS_67Cu = CS_67Zn_1*abund_67Zn + CS_67Zn_6*abund_67Zn + CS_68Zn_5*abund_68Zn + CS_70Zn_4*abund_70Zn
        #E_67Cu = E_67Zn_1

        # to make the plot smoother
        tck = interpolate.splrep(E_Zn, CS_Zn, s=0) #interpalation - give it raw data of energy array and cross section array that we have calculated
        E = np.linspace(0, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)
        #plt.plot(E, CS, color='c', label='EMPIRE')


        #print(len(E))
        #print(len(CS))
        print('...........................')
        #plt.plot(E, xaxis, 'y--', label='TENDL')
        plt.plot(E, CS, 'y', label='EMPIRE-3.2.3')
        #plt.plot(E_67Cu, CS_67Cu, 'y--', label='EMPIRE')

        return E, CS, tck

#--------------------- Al ------------------------

    if foil == 'Al':
    #A = ['191', '193'] # stable iridium isotopes
        abund_27Al = 1.00
    #file_ending =
    #endre f_191I til de stabile isotopene i zink osv
    #v_64Zn = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-0' + filename_CS + '_coh' + '.txt'
    # 64Zn
        u_27Al_empire = path + '/../EMPIRECOH2/' + foil + '/27Al/' +  Z + '-' + filename_CS + '_empire' + '.txt'


        if os.path.isfile(u_27Al_empire): # Na-24
            #print("Ir 191 file: ",f_191Ir)
            print("u_27Al_Na-24 exists")
            CS_27Al = np.genfromtxt(u_27Al_empire, delimiter='\t', usecols=[1])
            E_27Al = np.genfromtxt(u_27Al_empire, delimiter='\t', usecols=[0])
        else:
            print("27Al_Na-24 file does not exist")
            CS_27Al = 0
            E_27Al =  0


        CS_Al = CS_27Al*abund_27Al

        tck = interpolate.splrep(E_27Al, CS_Al, s=0) #interpalation - give it raw data of energy array and cross section array that we have calculated
        E = np.linspace(0, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)
        plt.plot(E, CS, 'y', label='EMPIRE-3.2.3')

        return E, CS, tck

#------------------------- In -----------------------------

    if foil == 'In':
    #A = ['191', '193'] # stable iridium isotopes
        abund_113In = 0.429 ; abund_115In = 0.9571
    #file_ending =
    #endre f_191I til de stabile isotopene i zink osv
    #v_64Zn = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-0' + filename_CS + '_coh' + '.txt'
        # 113In
        u_113In_empire = path + '/../EMPIRECOH2/' + foil + '/113In/' +  Z + '-' + filename_CS + '_empire' + '.txt'
        # 115In
        u_115In_empire = path + '/../EMPIRECOH2/' + foil + '/115In/' +  Z + '-' + filename_CS + '_empire' + '.txt'



########################## 113In #############################

        if os.path.isfile(u_113In_empire):
            #print("Ir 191 file: ",f_191Ir)
            #print("u_113In exists")
            CS_113In = np.genfromtxt(u_113In_empire, delimiter='\t', usecols=[1])
            E_In = np.genfromtxt(u_113In_empire, delimiter='\t', usecols=[0])
        else:
            #print("113In file does not exist")
            CS_113In = 0
            E_113In =  0


########################## 115In #############################

        if os.path.isfile(u_115In_empire):
            #print("Ir 191 file: ",f_191Ir)
            #print("u_115In exists")
            CS_115In = np.genfromtxt(u_115In_empire, delimiter='\t', usecols=[1])
            E_In = np.genfromtxt(u_115In_empire, delimiter='\t', usecols=[0])
        else:
            #print("115In_In-112 file does not exist")
            CS_115In = 0
            E_115In =  0



        CS_In = CS_113In*abund_113In + CS_115In*abund_115In

        tck = interpolate.splrep(E_In, CS_In, s=0) #interpalation - give it raw data of energy array and cross section array that we have calculated
        E = np.linspace(0, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)
        plt.plot(E, CS, 'y', label='EMPIRE-3.2.3')

        return E, CS, tck


#------------------------- Y -----------------------------

    if foil == 'Y':
    #A = ['191', '193'] # stable iridium isotopes
        abund_89Y = 1.00 ;
    #file_ending =
    #endre f_191I til de stabile isotopene i zink osv
    #v_64Zn = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-0' + filename_CS + '_coh' + '.txt'
        # 89Y
        u_89Y_empire = path + '/../EMPIRECOH2/' + foil + '/89Y/' +  Z + '-' + filename_CS + '_empire' + '.txt'



        if os.path.isfile(u_89Y_empire): # Y-87
            #print("Ir 191 file: ",f_191Ir)
            #print("u_89Y exists")
            CS_89Y = np.genfromtxt(u_89Y_empire, delimiter='	', usecols=[1])
            E_Y = np.genfromtxt(u_89Y_empire, delimiter='	', usecols=[0])
        else:
            #print("89Y file does not exist")
            CS_89Y = 0
            E_89Y =  0


        CS_Y = CS_89Y*abund_89Y

        tck = interpolate.splrep(E_Y, CS_Y, s=0) #interpalation - give it raw data of energy array and cross section array that we have calculated
        E = np.linspace(0, 40, 1000)
        CS = interpolate.splev(E, tck, der=0)
        plt.plot(E, CS, 'y', label='EMPIRE-3.2.3')

        return E, CS, tck

#------------------------- Zr -----------------------------


    if foil == 'Zr':
    #A = ['191', '193'] # stable iridium isotopes
        abund_90Zr = 0.5145 ; abund_91Zr = 0.1122 ; abund_92Zr = 0.1715 ; abund_94Zr = 0.1738 ; abund_96Zr = 0.280
    #file_ending =
    #endre f_191I til de stabile isotopene i zink osv
    #v_64Zn = path + '/../EMPIRECOH2/' + foil + '/64Zn/' +  Z + '-0' + filename_CS + '_coh' + '.txt'

######################### 90 Zr ############################
        u_90Zr_empire = path + '/../EMPIRECOH2/' + foil + '/90Zr/' +  Z + '-' + filename_CS + '_empire' + '.txt'
        u_91Zr_empire = path + '/../EMPIRECOH2/' + foil + '/91Zr/' +  Z + '-' + filename_CS + '_empire' + '.txt'
        u_92Zr_empire = path + '/../EMPIRECOH2/' + foil + '/92Zr/' +  Z + '-' + filename_CS + '_empire' + '.txt'
        u_94Zr_empire = path + '/../EMPIRECOH2/' + foil + '/94Zr/' +  Z + '-' + filename_CS + '_empire' + '.txt'
        u_96Zr_empire = path + '/../EMPIRECOH2/' + foil + '/96Zr/' +  Z + '-' + filename_CS + '_empire' + '.txt'

        E_Zr = 0.0

        if os.path.isfile(u_90Zr_empire):
            #print("Ir 191 file: ",f_191Ir)
            #print("u_90Zr exists")
            CS_90Zr = np.genfromtxt(u_90Zr_empire, delimiter='	', usecols=[1])
            E_Zr = np.genfromtxt(u_90Zr_empire, delimiter='	', usecols=[0])
        else:
            #print("90Zr file does not exist")
            CS_90Zr = 0
            E_90Zr =  0


######################### 91 Zr ############################


        if os.path.isfile(u_91Zr_empire):
            #print("Ir 191 file: ",f_191Ir)
            #print("u_91Zr exists")
            CS_91Zr = np.genfromtxt(u_91Zr_empire, delimiter='	', usecols=[1])
            E_Zr = np.genfromtxt(u_91Zr_empire, delimiter='	', usecols=[0])
        else:
            #print("91Zr file does not exist")
            CS_91Zr = 0
            E_91Z1 =  0

######################### 92 Zr ############################


        if os.path.isfile(u_92Zr_empire): # Y-90
            #print("Ir 191 file: ",f_191Ir)
            #print("u_92Zr exists")
            CS_92Zr = np.genfromtxt(u_92Zr_empire, delimiter='	', usecols=[1])
            E_Zr = np.genfromtxt(u_92Zr_empire, delimiter='	', usecols=[0])
        else:
            #print("92Zr file does not exist")
            CS_92Zr = 0
            E_92Zr =  0


######################### 94 Zr ############################


        if os.path.isfile(u_94Zr_empire):
            #print("Ir 191 file: ",f_191Ir)
            #print("u_94Zr exists")
            CS_94Zr = np.genfromtxt(u_94Zr_empire, delimiter='	', usecols=[1])
            E_Zr = np.genfromtxt(u_94Zr_empire, delimiter='	', usecols=[0])
        else:
            #print("94Zr file does not exist")
            CS_94Zr = 0
            E_94Zr =  0



######################### 96 Zr ############################


        if os.path.isfile(u_96Zr_empire):
            #print("Ir 191 file: ",f_191Ir)
            #print("u_96Zr exists")
            CS_96Zr = np.genfromtxt(u_96Zr_empire, delimiter='	', usecols=[1])
            E_Zr = np.genfromtxt(u_96Zr_empire, delimiter='	', usecols=[0])
        else:
            #print("96Zr file does not exist")
            CS_96Zr = 0
            E_96Zr =  0


        CS_Zr = CS_90Zr*abund_90Zr + CS_91Zr*abund_91Zr + CS_92Zr*abund_92Zr + CS_94Zr*abund_94Zr + CS_96Zr*abund_96Zr
        print('length E_ZR :', len(E_Zr))
        print('CS_Zr :', CS_Zr)

        if len(E_Zr) > 1:
            tck = interpolate.splrep(E_Zr, CS_Zr, s=0) #interpalation - give it raw data of energy array and cross section array that we have calculated
            E = np.linspace(0, 40, 1000)
            CS = interpolate.splev(E, tck, der=0)
            plt.plot(E, CS, 'y', label='EMPIRE-3.2.3')
            return E, CS, tck
        else:
            tck = interpolate.splrep(np.linspace(0, 40), np.zeros(len(np.linspace(0,40))) , s=0)
            return E_Zr, CS_Zr, tck


def IRDFF(foil, filename_CS_me):
    print('------------------- IRDFF--------------')
    if foil == 'Zr':
        irdff = path + '/../IRDFF/' + '/Zr/' + filename_CS_me + '.txt'

    elif foil == 'Y':
        irdff  = path + '/../IRDFF/' + '/Y/' + filename_CS_me + '.txt'

    elif foil == 'Al':
        irdff  = path + '/../IRDFF/' + '/Al/' + filename_CS_me + '.txt'

    print(irdff)
    if os.path.isfile(irdff):

        with open(irdff) as f:
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
                #dE.append(float(string[1]))
                CS.append(float(string[1])*1e3) # in mb
                dCS.append(float(string[2])*1e3) # in mb
                #author.append(string[5]) #index 4 is equal to #

            #print('EXFOR CS: ', CS)
            #plt.plot(E, CS)
                #plt.show()
                #plt.errorbar(E[ind], CS[ind], marker='.', linewidth=0.001, xerr=dE[ind], yerr=dCS[ind], elinewidth=0.5, capthick=0.5, capsize=3.0, label=author[ind] )
            plt.plot(E, CS, label='Recommended CS (IRDFF-II)')
            plt.fill_between(E, np.array(CS)+np.array(dCS), np.array(CS)-np.array(dCS), color='blue', alpha=0.1)
            #plt.errorbar(E, CS, marker='.', markersize=4, linewidth=0.001, xerr=dE, yerr=dCS, elinewidth=0.5, capthick=0.5, capsize=3.0, label=author[0] )
            #plt.legend()
            #plt.show()

            return E, CS, dCS



def mydata(filename_CS, foil, plot_flag):
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
        if plot_flag == True:

            #plt.plot(CS_energi, xaxis, 'g.', label='THIS WORK')
            #plt.plot(CS_energi, CS_zn, 'k.', label='This Work')
            #plt.errorbar(CS_energi, CS_zn, marker='P', color='darkred',linewidth=0.0001, yerr=dCS, elinewidth=1.0, capthick=1.0, capsize=3.0, label='This Work')
            plt.errorbar(CS_energi, zero_to_nan(CS_zn), color='k', marker='.', markersize=8 ,linewidth=0.001, xerr=dE, yerr=dCS, elinewidth=0.5, capthick=0.5, capsize=3.0, label='This Work' )
            #plt.errorbar(CS_energi[1], CS_zn[1], marker='.', linewidth=0.001, xerr=dE[1], yerr=dCS[1], elinewidth=0.5, capthick=0.5, capsize=3.0, label='This Work' )
            #plt.errorbar(CS_energi, CS_zn, marker='.', linewidth=0.001, xerr=dE, yerr=dCS, elinewidth=0.5, capthick=0.5, capsize=3.0 )

    #plt.show()

    return CS_energi, CS_zn, dCS, dE


def zero_to_nan(values):
    """Replace every 0 with 'nan' and return a copy."""
    return [float('nan') if x==0 else x for x in values]



def Cross_section(foil, A, Z, reaction, filename_CS_me, filename_CS, y_max=None, file_ending='.tot'):

    with open('../../Jon code/meulders_33MeV.csv') as f:
    	meulders_33MeV = np.array([i.split(',') for i in f.read().split('\n')[:-1]], dtype=np.float64)

    with open('../../Jon code/meulders_16MeV.csv') as f:
    	meulders_16MeV = np.array([i.split(',') for i in f.read().split('\n')[:-1]], dtype=np.float64)

    E_16 = meulders_16MeV[:,0]
    E_33 = meulders_33MeV[:,0]
    dPhidE_16 = meulders_16MeV[:,1]
    dPhidE_33 = meulders_33MeV[:,1]


    E_ALICE, CS_ALICE, tck_ALICE = ALICE(foil, A, Z)
    if len(A) == 2:
       A = '0' + A
    if len(Z) == 2:
       Z = '0' + Z


    E_TALYS, CS_TALYS, tck_TALYS = TALYS(foil, A, Z)
    EXFOR(foil, reaction)
    #no EXFOR for :
    # natZn(n,x)67Cu
    E_TENDL, CS_TENDL, tck_TENDL = Tendl(foil, A, Z)
    #print('Foil, A, Z, filename :', foil, A, Z, filename_CS_me)
    E_CoH, CS_CoH, tck_CoH = CoH(foil, A, Z, filename_CS_me)
    E_EMPIRE, CS_EMPIRE, tck_EMPIRE = EMPIRE(foil, A, Z, filename_CS)

    if filename_CS_me == '89Zr' or filename_CS_me == '24Na' or filename_CS_me == '88Y':
        E_IRDFF, CS_IRDFF, dEC_IRDFF = IRDFF(foil, filename_CS_me)

# integral for each of the codes

    # plt.figure(0)
    # plt.plot(meulders_33MeV[:,0],meulders_33MeV[:,1])
    # plt.show()
    #E_16_average = np.trapz(res_list_16, E_16) / np.trapz(dPhidE_16, E_16)
    flux_avg_CS_ALICE = [np.trapz(interpolate.splev(E_16, tck_ALICE, der=0)*dPhidE_16, E_16) / np.trapz(dPhidE_16, E_16), np.trapz(interpolate.splev(E_33, tck_ALICE, der=0)*dPhidE_33, E_33) / np.trapz(dPhidE_33, E_33)]
    #print('flux_avg_CS_ALICE : ', flux_avg_CS_ALICE)
    flux_avg_CS_TALYS = [np.trapz(interpolate.splev(E_16, tck_TALYS, der=0)*dPhidE_16, E_16) / np.trapz(dPhidE_16, E_16), np.trapz(interpolate.splev(E_33, tck_TALYS, der=0)*dPhidE_33, E_33) / np.trapz(dPhidE_33, E_33)]
    flux_avg_CS_TENDL = [np.trapz(interpolate.splev(E_16, tck_TENDL, der=0)*dPhidE_16, E_16) / np.trapz(dPhidE_16, E_16), np.trapz(interpolate.splev(E_33, tck_TENDL, der=0)*dPhidE_33, E_33) / np.trapz(dPhidE_33, E_33)]
    flux_avg_CS_CoH = [np.trapz(interpolate.splev(E_16, tck_CoH, der=0)*dPhidE_16, E_16) / np.trapz(dPhidE_16, E_16), np.trapz(interpolate.splev(E_33, tck_CoH, der=0)*dPhidE_33, E_33) / np.trapz(dPhidE_33, E_33)]
    flux_avg_CS_EMPIRE = [np.trapz(interpolate.splev(E_16, tck_EMPIRE, der=0)*dPhidE_16, E_16) / np.trapz(dPhidE_16, E_16), np.trapz(interpolate.splev(E_33, tck_EMPIRE, der=0)*dPhidE_33, E_33) / np.trapz(dPhidE_33, E_33)]





    CS_energi_mydata, CS_zn_mydata, dCS_mydata, dE_mydata = mydata(filename_CS_me, foil, plot_flag=False)

    #plt.plot(CS_energi_mydata ,flux_avg_CS_ALICE, color='y')
    #plt.plot(CS_energi_mydata ,flux_avg_CS_TALYS, color='g')
    #plt.plot(CS_energi_mydata ,flux_avg_CS_TENDL, color='r')
    #plt.plot(CS_energi_mydata ,flux_avg_CS_CoH, color='c')
    #plt.plot(CS_energi_mydata ,flux_avg_CS_EMPIRE, 'y--')

    # edit from plot to errorbar - play with linewith and markersize to see the points better
    plt.errorbar(CS_energi_mydata, flux_avg_CS_ALICE, color='m', marker='.', linewidth=0.001, xerr=dE_mydata, elinewidth=0.5, capthick=0.5, capsize=3.0 )
    plt.errorbar(CS_energi_mydata, flux_avg_CS_TALYS, color='g', marker='.', linewidth=0.001, xerr=dE_mydata, elinewidth=0.5, capthick=0.5, capsize=3.0 )
    plt.errorbar(CS_energi_mydata, flux_avg_CS_TENDL, color='r', marker='.', linewidth=0.001, xerr=dE_mydata, elinewidth=0.5, capthick=0.5, capsize=3.0 )
    plt.errorbar(CS_energi_mydata, flux_avg_CS_CoH, color='c', marker='.', linewidth=0.001, xerr=dE_mydata, elinewidth=0.5, capthick=0.5, capsize=3.0 )
    plt.errorbar(CS_energi_mydata, flux_avg_CS_EMPIRE, color='y', marker='.', linewidth=0.001, xerr=dE_mydata, elinewidth=0.5, capthick=0.5, capsize=3.0 )

    CS_energi_mydata, CS_zn_mydata, dCS_mydata, dE_mydata = mydata(filename_CS_me, foil, plot_flag=True)

    plt.legend(loc='upper left')
    plt.title(reaction)
    plt.xlim(0,40)
    plt.ylim(ymin=0)
    plt.ylim(ymax=y_max)
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
# LINE 1507 - Cross_section


#Cross_section('foil', 'A', 'Z', 'EXFOR', 'mydata')
#Cross_section('foil', 'A', 'Z', 'reaction', 'filename_CS_me', filename_CS)

### ZINK ###

#Cross_section('Zn', '62', '30', 'Zn(n,x)62Zn', '62ZN', 'Zn-62', y_max=30)
#Cross_section('Zn', '63', '30', 'Zn(n,x)63Zn', '63ZN', 'Zn-63')
#Cross_section('Zn', '64', '29', 'Zn(n,x)64Cu', '64CU', 'Cu-64')
#Cross_section('Zn', '65', '28', 'Zn(n,x)65Ni', '65NI', 'Ni-65', y_max=8)
#Cross_section('Zn', '65', '30', 'Zn(n,x)65Zn', '65ZN', 'Zn-65')

# NOT REPORT 66Cu - not enough info due to half-life and the time we counted
####### Cross_section('Zn', '66', '29', 'Zn(n,x)66Cu', '66CU', 'Cu-66')
#Cross_section('Zn', '67', '29', 'Zn(n,x)67Cu', '67CU', 'Cu-67')
######## Cross_section('Zn', '66', '28', 'Zn(n,x)66Ni', '66NI', 'Ni-66') #NOT REPORT
#Cross_section('Zn', '69', '30', 'Zn(n,x)69mZn', '69ZNm', 'Zn-69M')


### Zirconium ###

#Cross_section('Zr', '90', '39', 'Zr(n,x)90mY', '90Ym', 'Y-90M')
#Cross_section('Zr', '91', '38', 'Zr(n,x)91Sr', '91SR', 'Sr-91', y_max=6)
#Cross_section('Zr', '91', '39', 'Zr(n,x)91mY', '91YM', 'Y-91M')
#Cross_section('Zr', '92', '39', 'Zr(n,x)92Y', '92Y', 'Y-92')
#Cross_section('Zr', '93', '39', 'Zr(n,x)93Y', '93Y', 'Y-93')  DO THIS LATER!!!!!!!!!!!!!!!!!!!!
#Cross_section('Zr', '95', '40', 'Zr(n,x)95Zr', '95Zr', 'Zr-95')
#Cross_section('Zr', '97', '40', 'Zr(n,x)97Zr', '97Zr', 'Zr-97', y_max=6)
#Cross_section('Zr', '89', '40', 'Zr(n,x)89Zr', '89Zr', 'Zr-89')
Cross_section('Zr', '93', '39', 'Zr(n,x)93Y', '95NB', 'Nb-95') #INGEN EXFOR ELLER TENDL
#Cross_section('Zr', '93', '39', 'Zr(n,x)93Y', '97NB', 'Nb-97') #INGEN EXFOR ELLER TENDL



### Indium ###

#Cross_section('In', '111', '49', 'In(n,x)111In', '111In', 'In-111')
#Cross_section('In', '112', '49', '113In(n,2n)112In', '112IN')
#Cross_section('In', '112', '49', 'In(n,x)112mIn', '112INm')
#Cross_section('In', '112', '49', 'In(n,x)112mIn', '112INm')


###  Aluminum ###

#Cross_section('Al', '24', '11', '27Al(n,x)24Na', '24Na', 'Na-24')

### Yttrium ###

#Cross_section('Y', '87', '39', 'Y(n,x)87Y', '87Y', 'Y-87')
#Cross_section('Y', '89', '39', '', '87Ym')
#Cross_section('Y', '89', '39', '89Y(n,g)90mY', '90Ym')
#Cross_section('Y', '89', '39', '', '87SRm')
