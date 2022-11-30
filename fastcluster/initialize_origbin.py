import numpy as np
import matplotlib.pyplot as plt
import os

import auxiliary as auxi
import constants as const
import classes as cl

##################INITIALIZE 1st GENERATION ORIGINAL BINARY#################
def initialize_origbin(BBH,SC,flags,fil,gf):

    N=len(BBH.m1)

    #BBH.idi=np.arange(0,len(m1),1)

    #########################################################################
    #############read m1, m2, tform, sma, ecc from input files###############
    #########################################################################

    
    m1o=[]
    m2o=[]
    tformo=[]
    smao=[]
    ecco=[]
    idio=[]
    
    #idii,tformi,m1ii,m2ii,smai,ecci=np.genfromtxt(fil,usecols=(0,1,2,3,4,5),unpack=True,dtype=float)
    m1ii,m2ii,smai,ecci,tSN1,tSN2=np.genfromtxt(fil,usecols=(3,4,7,8,19,20),skip_header=3,unpack=True,dtype=float)
    #check that second SN is in tSN2
    tf=np.where(tSN1>tSN2)
    tt=tf[0]
    tSN1[tt],tSN2[tt]=tSN2[tt],tSN1[tt]
    tformi=np.copy(tSN2)

    mBHaverage=np.mean(m1ii)
    max1g=max(max(m1ii),max(m2ii))
    print(len(m1ii),N)    
    a=np.random.randint(0,len(m1ii),N)
    m1o=m1ii[a]
    m2o=m2ii[a]
    #   tdelo=tdeli[a]
     
    b=np.where(m1o<m2o)
    print("b0=", b[0])
    bb=b[0]
    m1o[b[0]],m2o[b[0]]=m2o[b[0]],m1o[b[0]]
    tformo=tformi[a]
    smao=smai[a]
    ecco=ecci[a]


    #flags.hard=np.zeros(len(m1o),dtype=int)
    #for i in range(len(m1o)):
    flags.hard=auxi.check_hard(m1o,m2o,smao,SC.maverage,SC.sigma)
        
    aa=np.where(flags.hard==0)
    a=aa[0]
    print("Soft in cases= ", len(a))

    BBH.m1=np.copy(m1o)
    BBH.m2=np.copy(m2o)
    BBH.sma=np.copy(smao)
    BBH.ecc=np.copy(ecco)
    tform=np.copy(tformo)
    
    BBH.idi=np.arange(0,BBH.numBH,1)

    
    #BBH.numBH=len(a)

    #BBH.m1=np.copy(m1o[a])
    #BBH.m2=np.copy(m2o[a])
    ##print(len(BBH.m2))
    #BBH.sma=np.copy(smao[a])
    #BBH.ecc=np.copy(ecco[a])
    #tform=np.copy(tformo[a])

    #SC.numBH=len(a)
    #SC.SClifetime=SC.SClifetime[a]
    #SC.vesc=SC.vesc[a]
    #SC.maverage=SC.maverage[a]
    #SC.cdensity=SC.cdensity[a]
    #SC.feq=SC.feq[a]
    #SC.rho=SC.rho[a]
    #SC.rhocgs=SC.rhocgs[a]
    #SC.Mtot=SC.Mtot[a]
    #SC.sigma1D=SC.sigma1D[a]
    #SC.sigma=SC.sigma[a]
    #SC.mmax=SC.mmax[a]
    #SC.trh0=SC.trh0[a]
    #SC.rh=SC.rh[a]

    #flags=cl.merger_flags(len(a))
    
    #BBH.idi=np.arange(0,BBH.numBH,1)


    file_timescale=gf.out_met+"/timescales_1stgen_orig.txt"
    ft=open(file_timescale,"w")
    ft.write("t_form/Myr\n")
    for tt in range(len(tform)):
        ft.write(str(BBH.idi[tt])+" "+str(tform[tt])+"\n")
    ft.close()



    
    # CHECK IF SN kick > vesc from cluster - If this is the case
    # do not consider this black hole any further
    flags.flagSN=np.empty(len(BBH.m1),dtype=object)
    flags.flag_t3bb=np.empty(len(BBH.m1),dtype=object)
    countSN=0
    count3bb=0
    for i in range(BBH.numBH): 
        if(gf.SN_kick==True):
            vSN=auxi.SN_kick(BBH.m1[i]+BBH.m2[i]) #to be updated with cm velocity
            if(vSN>SC.vesc[i]):
                flags.flagSN[i]="ejected_by_SN"
                countSN+=1
            else:
                flags.flagSN[i]="non_ejected_by_SN"
        else:
                flags.flagSN[i]="non_ejected_by_SN"
                
        if(tform[i]>min(const.tHubble,SC.SClifetime[i])):
            flags.flag_t3bb[i]="t3bb_too_long"
            count3bb+=1
        else:
            flags.flag_t3bb[i]="t3bb_ok"

        BBH.chi1[i]=auxi.spin(gf.spinrms)
        BBH.chi2[i]=auxi.spin(gf.spinrms)
    print("Original BBHs ejected by supernova= ",countSN)
    print("Original BBHs with tform>tSC= ",count3bb)

    
    th1=np.random.rand(len(BBH.m1)) #spin tilt BH1 (0deg parallel with orb.ang.mom.)
    th2=np.random.rand(len(BBH.m1)) #spin tilt BH2 (0deg parallel with orb.ang.mom.)
    BBH.th1=auxi.set_theta(th1)
    BBH.th2=auxi.set_theta(th2)
    print("Print len hard= ", len(BBH.th1))

    
    if(gf.plotta==True):
        plt.hist(BBH.m1)
        plt.hist(BBH.m2)
        plt.legend(["m1","m2"])
        plt.xlabel("Mass [M$_\odot$]")
        plt.tight_layout()
        plt.show()

        binni=np.logspace(-3,4,num=50)
        #print(binni)
        plt.hist(tform,bins=binni,density=True)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("T$_\mathrm{form}$ [Myr]")
        plt.tight_layout()
        plt.show()

        plt.hist(BBH.sma)
        plt.xlabel("SMA [R$_\odot$]")
        plt.show()
    
        plt.hist(BBH.ecc)
        plt.xlabel("ecc")
        plt.show()


    
    return BBH, SC, flags, m1ii, mBHaverage, max1g, tform
