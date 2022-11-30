import numpy as np

import auxiliary as auxi


##################INITIALIZE N-TH GENERATION DYNAMICAL BINARY
def initialize_nthgen(BBH,SC,max1g,flags,gf): #initializes N-gen binary

    N = len(BBH.m1)
    
    
    for i in range(N):
        if(gf.SN_kick==True):
            vSN=1e20
            while(vSN>SC.vesc[i]):
                if(gf.ngng==True):
                    BBH.m2[i]=auxi.secondary(BBH.m1[i],0.9*2.*BBH.m2min)
                else:
                    BBH.m2[i]=auxi.secondary(BBH.m1[i],BBH.m2min)
                vSN=auxi.SN_kick(BBH.m2[i])
        else:
            if(gf.ngng==True):
                BBH.m2[i]=auxi.secondary(BBH.m1[i],0.9*2.*BBH.m2min)
            else:
                BBH.m2[i]=auxi.secondary(BBH.m1[i],BBH.m2min)
                

        if((gf.ngng==True) or (BBH.m2[i]>max1g)):
            pp=np.random.randint(0,N)
            BBH.chi2[i]=BBH.chi1[pp]
        else:
            BBH.chi2[i]=auxi.spin(gf.spinrms)




        #sma_ngen=auxi.semimajor(m1_ngen) #now totally arbitrary
        #auxi.check_hard(m1,m2,sma,maverage,sigma)
    
    BBH.sma=auxi.semimajor(BBH.m1,BBH.m2) #now totally arbitrary
    flags.hard=auxi.check_hard(BBH.m1,BBH.m2,BBH.sma,SC.maverage,SC.sigma)
    a=np.where(flags.hard==0)
    print("soft nth generation binaries= ",len(a[0]))
    
    BBH.ecc=auxi.ecc(BBH.m1) #now totally arbitrary

    th1_ngen=np.random.rand(N) #spin tilt BH1 (0deg parallel with orb.ang.mom.)
    th2_ngen=np.random.rand(N) #spin tilt BH2 (0deg parallel with orb.ang.mom.)

    BBH.th1=auxi.set_theta(th1_ngen)
    BBH.th2=auxi.set_theta(th2_ngen)
    
    return BBH 
