import numpy as np

import peters as pt
import peters_exchanges as ptexch

import auxiliary as auxi
import constants as const

###################CORE CALCULATION###################################
def make_generation_evol(BBH, SC, t3bb, vesc0, rh0, mBHaverage, m1ii, flags, gf, NC):
    #if dynamical binary and ejected by the supernova before formation do not consider it any further (to be updated when stellar evolution included)

    N=len(BBH.m1)


    
    tdel=np.zeros(N,float)
    tmerg=np.zeros(N,float)


    print(len(SC.maverage),len(BBH.m1),len(BBH.m2),len(SC.vesc))
    aej=auxi.a_ej(BBH.xi,SC.maverage,BBH.m1,BBH.m2,SC.vesc)  #semi-major axis for dynamical ejection
    aGW=auxi.a_GW(BBH.ecc,BBH.m1,BBH.m2,BBH.xi,SC.sigma,SC.rhocgs*const.c2hm_dens) #semi-major axis where GWs dominate
    print(aej[0]/const.Rsun)
    print(aGW[0]/const.Rsun)

    for i in range(len(BBH.m1)):

        #do not integrate soft binaries
        if(flags.hard[i]==0): 
            flags.flag[i]="soft_binary"
            flags.flag2[i]="no_merg"
            flags.flag3[i]="soft_binary"
            flags.flag_evap[i]="soft_binary"
            tmerg[i]=10.*const.tHubble
            continue
        
        #do not integrate if t3bb or tform is too long
        if(flags.flag_t3bb[i]=="t3bb_too_long"): 
            flags.flag[i]="t3bb_too_long"
            flags.flag2[i]="no_merg"
            flags.flag3[i]="t3bb_too_long"
            flags.flag_evap[i]="t3bb_too_long"
            tmerg[i]=10.*const.tHubble
            continue

        #make sure dynamical binary which does not form because of SN kick is not integrated
        if((gf.channel=="Dyn") and (flags.flagSN[i]=="ejected_by_SN")):
            flags.flag[i]="ejected_by_SN"
            flags.flag2[i]="no_merg"
            flags.flag3[i]="ejected_by_SN"
            flags.flag_evap[i]="ejected_by_SN"
            tmerg[i]=10.*const.tHubble
            continue

        if(aej[i]>aGW[i]):
            flags.flag[i]="ejected_3B" #ejected by 3body encounters
        else:
            flags.flag[i]="in_cluster" #non ejectedby 3body encounters
            #print(flags.flag)


        if(gf.exchanges==True):
            BBH, SC, tdel[i], flags = ptexch.peters_evolv(BBH, SC, i, t3bb[i], vesc0[i], rh0[i], mBHaverage, m1ii, aej[i], aGW[i], flags) #integrate GW + scattering
        else:
            BBH, SC, tdel[i], flags = pt.peters_evolv(BBH, SC, i, t3bb[i], vesc0[i], rh0[i], aej[i], aGW[i], flags, NC) #integrate GW + scattering
 

        tdel[i]=tdel[i]/(1e6*const.yr)
        tmerg[i]=t3bb[i]+tdel[i]


        #make sure original binary ejected by SN kick is not retained for next generation
        if((gf.channel=="Orig") and (flags.flagSN[i]=="ejected_by_SN")):
            flags.flag[i]="ejected_by_SN"
            flags.flag3[i]="ejected_by_SN"
            flags.flag_evap[i]="ejected_by_SN"
            
            if(flags.flag2[i]=="no_merg"): #flag2 assigned by peters
                tmerg[i]=10.*const.tHubble
                flags.flag3[i]=flags.flag2[i]
                
            elif(flags.flag2[i]=="merg"):
                BBH.mmerg[i]=auxi.final_mass(BBH.m1[i],BBH.m2[i],BBH.chi1[i],BBH.chi2[i])
                BBH.chimerg[i]=auxi.final_spin(BBH.chi1[i], BBH.chi2[i], BBH.m1[i], BBH.m2[i])
            else:
                print("error in orig ejected")

            continue


        
        #all binary systems not ejected by SN and with short t3bb/tform
        if((flags.flag_t3bb[i]=="t3bb_ok") and (flags.flagSN[i]=="non_ejected_by_SN")):
            if(flags.flag2[i]=="no_merg"): #flag2 assigned by peters
                tmerg[i]=10.*const.tHubble
                flags.flag3[i]=flags.flag2[i]
                
            elif(flags.flag2[i]=="merg"):
                BBH.mmerg[i]=auxi.final_mass(BBH.m1[i],BBH.m2[i],BBH.chi1[i],BBH.chi2[i])
                BBH.chimerg[i]=auxi.final_spin(BBH.chi1[i], BBH.chi2[i], BBH.m1[i], BBH.m2[i])
                if(flags.flag[i]=="ejected_3B"):
                    flags.flag3[i]=flags.flag[i] #ejected by 3body encounters
                elif(flags.flag[i]=="in_cluster"):
                    if((flags.flag_evap[i]=="cluster_died") or (flags.flag_evap[i]=="ejected_3B")):
                        flags.flag3[i]=flags.flag_evap[i]
                    elif(flags.flag_evap[i]=="in_cluster"):
                        BBH.vk[i]=auxi.relativistic_kick(BBH.m1[i],BBH.m2[i],BBH.chi1[i],BBH.chi2[i],BBH.th1[i],BBH.th2[i])
                        if(BBH.vk[i]<=SC.vesc[i]): 
                            flags.flag3[i]="in_cluster"
                        else:
                            flags.flag3[i]="ejected_GR" #ejected either by GR
                    else:
                        flags.flag3[i]="error"
                        #print(flags.flag3[i],flags.flag_evap[i])
                else:
                    flags.flag3[i]="error"
                    #print(flags.flag3[i],flags.flag[i])
            else:
                flags.flag3[i]="error"
                print(flags.flag3[i],flags.flag2[i])
        else:
            flags.flag3[i]="error"
            #print(flags.flag3[i],flags.flag_t3bb[i])

    SC.rhocgs=np.copy(SC.rho*const.msun/(const.parsec**3.))

    SC.sigma1D=SC.sigma/np.sqrt(3.)

    return BBH, SC, tdel, tmerg, flags





def make_generation(BBH, SC, t3bb, mBHaverage, m1ii, flags, gf):

    N=len(BBH.m1)
    tdel=np.zeros(N,float)
    tmerg=np.zeros(N,float)



    aej=auxi.a_ej(BBH.xi,SC.maverage,BBH.m1,BBH.m2,SC.vesc) #2.*xi*maverage/(m1+m2)**3. * G * m1 * m2/vesc/vesc #semi-major axis for dynamical ejection
    aGW=auxi.a_GW(BBH.ecc,BBH.m1,BBH.m2,BBH.xi,SC.sigma,SC.rhocgs*const.c2hm_dens)
    print(aej[0]/const.Rsun)
    print(aGW[0]/const.Rsun)

    for i in range(len(BBH.m1)):

        #do not integrate soft binaries
        if(flags.hard[i]==0): 
            flags.flag[i]="soft_binary"
            flags.flag2[i]="no_merg"
            flags.flag3[i]="soft_binary"
            flags.flag_evap[i]="soft_binary"
            tmerg[i]=10.*const.tHubble
            continue
        
        #do not integrate if t3bb or tform is too long
        if(flags.flag_t3bb[i]=="t3bb_too_long"): 
            flags.flag[i]="t3bb_too_long"
            flags.flag2[i]="no_merg"
            flags.flag3[i]="t3bb_too_long"
            flags.flag_evap[i]="t3bb_too_long"
            tmerg[i]=10.*const.tHubble
            continue
        
        #make sure dynamical binary which does not form because of SN kick is not integrated
        if((gf.channel=="Dyn") and (flags.flagSN[i]=="ejected_by_SN")):
            flags.flag[i]="ejected_by_SN"
            flags.flag2[i]="no_merg"
            flags.flag3[i]="ejected_by_SN"
            flags.flag_evap[i]="ejected_by_SN"
            tmerg[i]=10.*const.tHubble
            #print("found ",flags.flagSN[i]," ",flags.flag3[i],i)
            continue
       
        if(aej[i]>aGW[i]):
            flags.flag[i]="ejected_3B"
        else:
            flags.flag[i]="in_cluster"
            #print(flags.flag)

        if(gf.exchanges==True):
            BBH, tdel[i], flags = ptexch.peters(BBH, SC, i, t3bb[i], mBHaverage, m1ii, aej[i], aGW[i], flags) # integrate GW+hardening
        else:
            BBH, tdel[i], flags = pt.peters(BBH, SC, i, t3bb[i], aej[i], aGW[i], flags) #m1[i],m2[i],sma[i],ecc[i],xi,ki,rhocgs[i],sigma[i],t3bb[i],SClifetime[i],aej[i],aGW[i]) # integrate GW+hardening


        tdel[i]/=(1e6*const.yr)
        tmerg[i]=t3bb[i]+tdel[i]


        #make sure original binary ejected by SN kick is not retained for next generation
        if((gf.channel=="Orig") and (flags.flagSN[i]=="ejected_by_SN")):
            flags.flag[i]="ejected_by_SN"
            flags.flag3[i]="ejected_by_SN"
            flags.flag_evap[i]="ejected_by_SN"
            
            if(flags.flag2[i]=="no_merg"): #flag2 assigned by peters
                tmerg[i]=10.*const.tHubble
                flags.flag3[i]=flags.flag2[i]
                
            elif(flags.flag2[i]=="merg"):
                BBH.mmerg[i]=auxi.final_mass(BBH.m1[i],BBH.m2[i],BBH.chi1[i],BBH.chi2[i])
                BBH.chimerg[i]=auxi.final_spin(BBH.chi1[i], BBH.chi2[i], BBH.m1[i], BBH.m2[i])
            else:
                print("error in orig ejected")
                
            continue


        
        #all binary systems not ejected by SN and with short t3bb/tform
        if((flags.flag_t3bb[i]=="t3bb_ok") and (flags.flagSN[i]=="non_ejected_by_SN")):
            if(flags.flag2[i]=="no_merg"): #flag2 assigned by peters
                tmerg[i]=10.*const.tHubble
                flags.flag3[i]=flags.flag2[i]
            elif(flags.flag2[i]=="merg"):
                BBH.mmerg[i]=auxi.final_mass(BBH.m1[i],BBH.m2[i],BBH.chi1[i],BBH.chi2[i])
                BBH.chimerg[i]=auxi.final_spin(BBH.chi1[i], BBH.chi2[i], BBH.m1[i], BBH.m2[i])
                if(flags.flag[i]=="ejected_3B"):
                    flags.flag3[i]=flags.flag[i] #ejected by 3body encounters
                elif(flags.flag[i]=="in_cluster"):
                    if((flags.flag_evap[i]=="cluster_died") or (flags.flag_evap[i]=="ejected_3B")):
                        flags.flag3[i]=flags.flag_evap[i]
                    elif(flags.flag_evap[i]=="in_cluster"):
                        BBH.vk[i]=auxi.relativistic_kick(BBH.m1[i],BBH.m2[i],BBH.chi1[i],BBH.chi2[i],BBH.th1[i],BBH.th2[i])
                        if(BBH.vk[i]<=SC.vesc[i]): 
                            flags.flag3[i]="in_cluster"
                        else:
                            flags.flag3[i]="ejected_GR" #ejected either by GR
                    else:
                        flags.flag3[i]="error"
                        #print(flags.flag3[i],flags.flag_evap[i])
                else:
                    flags.flag3[i]="error"
                    print(flags.flag3[i],flags.flag[i])
            else:
                flags.flag3[i]="error"
                #print(flags.flag3[i],flags.flag2[i])
        else:
            flags.flag3[i]="error"
            #print(flags.flag3[i],flags.flag_t3bb[i])


    return BBH, tdel, tmerg, flags


