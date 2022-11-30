import numpy as np
import matplotlib.pyplot as plt


import auxiliary as auxi
import constants as const


##################INITIALIZE 1st GENERATION DYNAMICAL BINARY#################
def initialize_dynbin(BBH, SC, flags, fil, gf, NC):
    N=len(BBH.m1)
    
    BBH.idi=np.arange(0,N,1)

    #########################################################################
    #############read primary mass and tform from input files################
    #########################################################################
    

    m1iii,tformii=np.genfromtxt(fil,usecols=(2,3),skip_header=3,unpack=True)
    if(gf.Z_list=="0.02"): #this if should be removed when bug in MOBSE fixed
        wrong=np.where((m1iii<=30.0) & (tformii<100.)) #this if should be removed when bug in MOBSE fixed
        ok=wrong[0]
        m1ii=np.copy(m1iii[ok])
        tformi=np.copy(tformii[ok])
    if(gf.Z_list!="0.02"): #this if should be removed when bug in MOBSE fixed
        wrong=np.where((m1iii<=80.0) & (tformii<100.)) #this if should be removed when bug in MOBSE fixed
        ok=wrong[0]
        m1ii=np.copy(m1iii[ok])
        tformi=np.copy(tformii[ok])
    print("Removed ",len(m1iii)-len(m1ii)," black holes")

    mBHaverage=np.mean(m1ii) #average BH mass of first generation
    max1g=max(m1ii) #maximum mass of first generation
    
    #print(len(m1ii),N)    


    a=np.random.randint(0,len(m1ii),N)
    BBH.m1=m1ii[a]
    tform=tformi[a]
    print(min(tform))
    #print(m1ii[a[100]],a[100],BBH.m1[100])

            
    #print("len of BBH is ",len(BBH.m1),len(m1ii))


    print("Maximum mass for 1g black holes is ",max1g)

    ###################################################################
    ###################CALCULATES t_DF, t3bb, t12######################
    ###################################################################
    t_DF=auxi.t_dynfric(SC.sigma,BBH.m1,SC.rhocgs)
    #t_DF=auxi.timescales(NC)[2]
    #dyn fric timescale in s
    t_DF=t_DF/const.yr/1e6 #dyn fric timescale in Myr
    #t_DF=auxi.t_dynfric2(m1,maverage,Mtot,rh) 

    t3bb=auxi.t_3bb(SC.cdensity,SC.feq,SC.sigma1D,BBH.m1) #125.0*(1e6/cdensity)**2  * (feq * sigma1D/30.e5)**9 * (20./m1)**5 #time for three-body capture in Myr, Lee (1995)
    hard=const.G * SC.maverage * const.msun/(SC.sigma* SC.sigma)
    t12init=auxi.t_12init(BBH.m1, SC.maverage, hard, SC.sigma, SC.cdensity, SC.binfrac)

    file_timescale=gf.out_met+"/timescales_1stgen.txt"
    ft=open(file_timescale,"w")
    ft.write("ID t_DF/Myr t_3bb/Myr t_12/Myr t_form/Myr\n")
    for tt in range(N):
        ft.write(str(tt)+" "+str(t_DF)+" "+str(t3bb[tt])+" "+str(t12init[tt])+" "+str(tform[tt])+"\n")
    ft.close()

    ##############################################################################################
    ##############calculates final tdyn (saved in t3bb) and other binary properties################
    ###############################################################################################
    count3bb=int(0)
    countSN=int(0)
    
    for i in range(N):
        t3bb[i]=min(t3bb[i],t12init[i])
        #print(t3bb[i],t_DF[i])
        t3bb[i]=t3bb[i]+t_DF[i]
        t3bb[i]=max(t3bb[i],tform[i]) #avoid that dynamical binary black hole forms before black hole has formed
        if(t3bb[i]>min(SC.SClifetime[i], const.tHubble)):
            flags.flag_t3bb[i]="t3bb_too_long"
            count3bb+=1
            #print(flags.flag_t3bb[i],t3bb[i],m1[i])
        else:
            flags.flag_t3bb[i]="t3bb_ok"

            
        BBH.m2[i]=auxi.secondary(BBH.m1[i],BBH.m2min)

        
        if(gf.SN_kick==True):
            # CHECK IF SN kick > vesc from cluster - If this is the case
            # do not consider this black hole any further
            vSN1=auxi.SN_kick(BBH.m1[i])
            vSN2=auxi.SN_kick(BBH.m2[i])
            #print("SN kick", vSN1,vSN2)
            if((vSN1>SC.vesc[i]) or (vSN2>SC.vesc[i])):
                flags.flagSN[i]="ejected_by_SN"
                countSN+=1
            else:
                flags.flagSN[i]="non_ejected_by_SN"
        else:
            flags.flagSN[i]="non_ejected_by_SN"
            




        BBH.chi1[i]=auxi.spin(gf.spinrms) #np.random.rand())
        BBH.chi2[i]=auxi.spin(gf.spinrms) #np.random.rand())



    BBH.sma=auxi.semimajor(BBH.m1,BBH.m2) #now totally arbitrary
    flags.hard=auxi.check_hard(BBH.m1,BBH.m2,BBH.sma,SC.maverage,SC.sigma)
    pr=np.where(flags.hard==0)
    print("The soft binaries are ",len(pr[0]))
    BBH.ecc=auxi.ecc(BBH.m1) #now totally arbitrary

    th1=np.random.rand(N) #spin tilt BH1 (0deg parallel with orb.ang.mom.)
    th2=np.random.rand(N) #spin tilt BH2 (0deg parallel with orb.ang.mom.)
    BBH.th1=auxi.set_theta(th1)
    BBH.th2=auxi.set_theta(th2)

    

    print("t3bb too long in ",count3bb," cases")
    print("ejected by SN kick in ",countSN," cases")


    return BBH,SC,flags,m1ii,mBHaverage, max1g, tform, t_DF, t3bb, t12init
