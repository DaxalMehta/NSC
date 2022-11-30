#############################################################################
########## FASTCLUSTER CODE V1.0 by Michela Mapelli #########################
########## CALCULATE MERGER FAMILY OF DYNAMICAL and ORIGINAL BINARIES #######
########## Release: April 25 2021 ###########################################
#############################################################################

import numpy as np
import matplotlib.pyplot as plt
import os

import classes as cl

import initialize_origbin as initorig
import initialize_dynbin as initdyn
import initialize_nthgen as initngen
import make_generation as makegen
import auxiliary as auxi
import constants as const


#### initialize the general flags
gf = cl.general_flags()
print(gf.out_top)


#### creates output path (if needed)
if not os.path.exists(gf.out_top):
    os.mkdir(gf.out_top)
    
if not os.path.exists(gf.out_dir):
    os.mkdir(gf.out_dir)

if not os.path.exists(gf.out_met):
    os.mkdir(gf.out_met)
   

    

#### initialize BBH properties
numBH=int(gf.Nref) 
BBH = cl.BBH_properties(numBH)
    

#### initialize SC properties
time=1000 # Initial time of simulation Myr
SC=cl.SC_properties(gf.flagNSC,int(numBH))
SC.init_SC(gf.flagNSC,int(numBH),time)
#print(SC.vesc)
vesc0=np.copy(SC.vesc) #save vesc0 for later

rh0=np.copy(SC.rh)     #save rh0 for later

#### initialize merger flags
flags = cl.merger_flags(len(BBH.m1))

### initialize nsc
NC=cl.gal_properties()
if(gf.channel=="Dyn"):
    fil=gf.in_dir+"/single_BH_"+str(gf.Z_list)+".txt"

    # initialize dynamical binary
    BBH,SC,flags,m1ii,mBHaverage, max1g, tform, t_DF, t3bb, t12init = initdyn.initialize_dynbin(BBH, SC, flags, fil, gf, NC)
        
elif(gf.channel=="Orig"):
    fil=gf.in_dir+"/data_BBHs_"+str(gf.Z_list)+".txt"

    # initialize original binary
    BBH,SC,flags, m1ii, mBHaverage, max1g, t3bb = initorig.initialize_origbin(BBH, SC, flags, fil, gf)
    #print(len(SC.SClifetime),len(BBH.m1))

elif(gf.channel!="Dyn" and gf.channel!="Orig"):
    print("something wrong, channel not found\n")
    exit()
        



##### print some diagnostics (if gf.Plotta==True)        
print("vesc logmean= ",10.**np.mean(np.log10(SC.vesc/1e5)))
if(gf.plotta==True):
    binss=np.logspace(np.log10(min(SC.vesc/1e5)),np.log10(max(SC.vesc/1e5)),num=100)
    plt.hist(SC.vesc/1e5,histtype="step", density=False, bins=binss)
    plt.xscale("log")
    plt.xlabel("V$_{esc}$ [km s$^{-1}$]")
    plt.tight_layout()
    plt.show()
    binss=np.logspace(np.log10(min(SC.Mtot)),np.log10(max(SC.Mtot)),num=100)
    plt.hist(SC.Mtot, histtype="step", density=False, bins=binss)
    plt.xscale("log")
    plt.xlabel("M$_{tot}$ [M$_{\odot}$]")
    plt.tight_layout()
    plt.show()
    binss=np.logspace(np.log10(min(SC.rho)),np.log10(max(SC.rho)),num=100)
    plt.hist(SC.rho,histtype="step", density=False, bins=binss)
    plt.xscale("log")
    plt.xlabel("$\\rho $ [M$_\odot$ pc$^{-3}$]")
    plt.tight_layout()
    plt.show()
    binss=np.logspace(np.log10(min(m1ii)),np.log10(max(m1ii)),num=100)
    plt.hist(m1ii,histtype="step", density=False, bins=binss)
    plt.xscale("log")
    plt.xlabel("$m_1$ [M$_\odot$]")
    plt.tight_layout()
    plt.show()
    binss=np.logspace(np.log10(min(t3bb)),np.log10(max(t3bb)),num=100)
    plt.hist(t3bb,histtype="step", density=False, bins=binss)
    plt.xscale("log")
    if(gf.channel=="Dyn"):
        plt.xlabel("$t_\mathrm{3bb}$ [Myr]")
    elif(gf.channel=="Orig"):
        plt.xlabel("$t_\mathrm{form}$ [Myr]")
    plt.tight_layout()
    plt.show()


    plt.hist(BBH.m1)
    plt.hist(BBH.m2)
    plt.legend(["m1","m2"])
    plt.xlabel("Mass [M$_\odot$]")
    plt.tight_layout()
    plt.show()
    
    plt.hist(BBH.sma)
    plt.xlabel("SMA [R$_\odot$]")
    plt.show()
    
    plt.hist(BBH.ecc)
    plt.xlabel("ecc")
    plt.show()



#############calculate and print first generation ##########################
Ngen = np.ones(len(BBH.m1))
    
if(gf.evol_cluster==True):
    BBH, SC, tdel, tmerg, flags = makegen.make_generation_evol(BBH, SC, t3bb, vesc0, rh0, mBHaverage, m1ii, flags, gf, NC)

else:
    BBH, tdel, tmerg, flags = makegen.make_generation(BBH, SC, t3bb, mBHaverage, m1ii, flags, gf)

    
out_file1=gf.out_met+"/first_gen.txt"
f=open(out_file1,"w")
f.write("c0:identifier c1:M1/Msun c2:M2/Msun c3:chi1 c4:chi2 c5:theta1 c6:theta2 c7:SMA(Rsun) c8:ecc c9:tDF+min(t12,t3bb)/Myr c10:SMAfin(cm) c11:eccfin c12:tpeters/Myr c13:(tDF+t3bb+tpeters)/Myr c14:vkick/kms c15:mrem/Msun c16:arem c17:vesc/kms c18:flag1 c19:flag2 c20:flag3 c21:flagSN c22:flag_exch c23:flag_t3bb c24:flag_evap c25:Mtot/Msun c26:ecc(10Hz) c27:Ngen\n")

if(gf.evol_cluster==True):
    for i in range(len(BBH.m1)):
        f.write(str(BBH.idi[i])+" "+str(BBH.m1[i])+" "+str(BBH.m2[i])+" "+str(BBH.chi1[i])+" "+str(BBH.chi2[i])+" "+str(BBH.th1[i])+" "+str(BBH.th2[i])+" "+str(BBH.sma[i])+" "+str(BBH.ecc[i])+" "+str(t3bb[i])+" "+str(BBH.sma2[i])+" "+str(BBH.ecc2[i])+" "+str(tdel[i])+" "+str(tmerg[i])+" "+str(BBH.vk[i]/1e5)+" "+str(BBH.mmerg[i])+" "+str(BBH.chimerg[i])+" "+str(vesc0[i]/1e5)+" "+str(flags.flag[i])+" "+str(flags.flag2[i])+" "+str(flags.flag3[i])+" "+str(flags.flagSN[i])+" "+str(flags.flag_exch[i])+" "+str(flags.flag_t3bb[i])+" "+str(flags.flag_evap[i])+" "+str(SC.Mtot[i])+" "+str(BBH.ecc10[i])+" "+str(Ngen[i])+"\n")
else:
    for i in range(len(BBH.m1)):
        f.write(str(BBH.idi[i])+" "+str(BBH.m1[i])+" "+str(BBH.m2[i])+" "+str(BBH.chi1[i])+" "+str(BBH.chi2[i])+" "+str(BBH.th1[i])+" "+str(BBH.th2[i])+" "+str(BBH.sma[i])+" "+str(BBH.ecc[i])+" "+str(t3bb[i])+" "+str(BBH.sma2[i])+" "+str(BBH.ecc2[i])+" "+str(tdel[i])+" "+str(tmerg[i])+" "+str(BBH.vk[i]/1e5)+" "+str(BBH.mmerg[i])+" "+str(BBH.chimerg[i])+" "+str(SC.vesc[i]/1e5)+" "+str(flags.flag[i])+" "+str(flags.flag2[i])+" "+str(flags.flag3[i])+" "+str(flags.flagSN[i])+" "+str(flags.flag_exch[i])+" "+str(flags.flag_t3bb[i])+" "+str(flags.flag_evap[i])+" "+str(SC.Mtot[i])+" "+str(BBH.ecc10[i])+" "+str(Ngen[i])+"\n")
f.close()

#### print some diagnostics (if Plotta == True)
if gf.plotta == True:
    plt.hist(BBH.vk/1e5)
    plt.xlabel("V$_k$ [km/s]")
    plt.show()

    a=np.where((flags.flag2=="merg") & (flags.flag_t3bb=="t3bb_ok")) #all merging systems, inside and outside cluster
    v=a[0]

    plt.hist(BBH.mmerg[v])
    plt.xlabel("M$_\mathrm{merg}$ [M$_\odot$]")
    plt.show()

    plt.hist(BBH.chimerg[v])
    plt.xlabel("a$_\mathrm{merg}$")
    plt.show()
    if gf.evol_cluster == True:
        plt.hist(vesc0/1e5)
        plt.hist(SC.vesc/1e5)
        plt.legend(["V$_{esc,0}$", "V$_{esc}$"])
        plt.xlabel("V$_{esc}$ [km/s]")
        plt.show()




####################################################################
#############CALCULATE AND PRINT N-TH GENERATIONS####################
#####################################################################

out_file2=gf.out_met+"/nth_generation.txt"
f2=open(out_file2,"w")
f2.write("c0:identifier c1:M1/Msun c2:M2/Msun c3:chi1 c4:chi2 c5:theta1 c6:theta2 c7:SMA(Rsun) c8:ecc c9:(tDF+min(t12,t3bb))/Myr c10:SMAfin(cm) c11:eccfin c12:tpeters/Myr c13:(ngen (tDF+t3bb+tpeters))/Myr c14:vkick/kms c15:mrem/Msun c16:arem c17:vesc/kms c18:flag1 c19:flag2 c20:flag3 c21:flagSN c22:flag_exch c23:flag_t3bb c24:flag_evap c25:Mtot/Msun c26:ecc(10Hz) c27:Ngen\n")

file_timescale=gf.out_met+"/timescales_nthgen.txt"
ft=open(file_timescale,"w")
ft.write("ID t_DF/Myr t_3bb/Myr t_12/Myr\n")

a=np.zeros(len(BBH.m1),int)
print(flags.flag3)
gen=2

while(len(a)>0):
    aa=np.where(flags.flag3=="in_cluster") # & (flags.flagSN=="non_ejected_by_SN")) #only mergers IN CLUSTER are considered for next generation
    a=np.copy(aa[0])
    print("Length of a", len(a))
    if(len(a)==0):
        print("last gen is ", gen-1)
        print(a)
        break
    else:
        m1_ngen=[]
        chi1_ngen=[]
        t3bb_ngen=[]
        idi_ngen=[]
        tmerg_ngen=[]
        Ngen=[]
        rhocgs_ngen=[]
        sigma_ngen=[]
        maverage_ngen=[]
        vesc_ngen=[]
        flagSN_ngen=[]
        flag_evap_ngen=[]
        flag_t3bb_ngen=[]
        Mtot_ngen=[]
        SClifetime_ngen=[]
        feq_ngen=[]
        binfrac_ngen=[]
        rho_ngen=[]
        rh_ngen=[]
        trh0_ngen=[]
        
        if(gf.evol_cluster==True):
            vesc0_ngen=[]
            rh0_ngen=[]



        for i in range(len(a)):
            #print(a[i])
            l=a[i]
            t_DF=auxi.t_dynfric(SC.sigma[l],BBH.mmerg[l],SC.rhocgs[l])
            #t_DF=auxi.timescales(NC)[2] 
            #dyn fric timescale in s

            t_DF=t_DF/const.yr/1e6 #dyn fric timescale in Myr

            #t_DF=auxi.t_dynfric2(mmerg[l],maverage[l],Mtot[l],rh[l]) 

            dens= const.c2hm_dens * SC.rho[l]/SC.maverage[l]
            sigma1=SC.sigma[l]/np.sqrt(3.)
            t3bb=auxi.t_3bb(dens,SC.feq[l],sigma1,BBH.mmerg[l]) #MM on april 2021 11 to change mmerg to average BH mass? note diff normal antonini/fragione
            
            hard=const.G * SC.maverage[l] * const.msun/(SC.sigma[l] * SC.sigma[l])
            t12init=auxi.t_12init(BBH.m1[l],SC.maverage[l],hard,SC.sigma[l],dens,SC.binfrac[l])
            #print("Times = ", t_DF,t3bb,t12init)
            ft.write(str(BBH.idi[l])+" "+str(t_DF)+" "+str(t3bb)+" "+str(t12init)+"\n")
            t3bb=min(t3bb,t12init)

            t3bb=t_DF+t3bb
            
            #print("TDF: ", t_DF)
            #print(tmerg[l],t3bb)
            if((tmerg[l]+t3bb)<min(SC.SClifetime[l],const.tHubble)):

                Ngen.append(gen)

                idi_ngen.append(BBH.idi[l])
                    
                m1_ngen.append(BBH.mmerg[l])
                chi1_ngen.append(BBH.chimerg[l])

                tmerg_ngen.append(tmerg[l])
                t3bb_ngen.append(t3bb)
                    
                flagSN_ngen.append(flags.flagSN[l])
                flag_evap_ngen.append(flags.flag_evap[l])
                flag_t3bb_ngen.append(flags.flag_t3bb[l])
                    
                rhocgs_ngen.append(SC.rhocgs[l])
                sigma_ngen.append(SC.sigma[l])
                maverage_ngen.append(SC.maverage[l])
                vesc_ngen.append(SC.vesc[l])


                Mtot_ngen.append(SC.Mtot[l])
                SClifetime_ngen.append(SC.SClifetime[l])
                feq_ngen.append(SC.feq[l])
                binfrac_ngen.append(SC.binfrac[l])
                rho_ngen.append(SC.rho[l])
                    
                trh0_ngen.append(SC.trh0[l])
                rh_ngen.append(SC.rh[l])
                
                if(gf.evol_cluster==True):
                    vesc0_ngen.append(vesc0[l])
                    rh0_ngen.append(rh0[l])



                        
        print("print len of m1", len(m1_ngen))

            
        BBH=cl.BBH_properties(len(m1_ngen))
            
        SC=cl.SC_properties(gf.flagNSC,len(m1_ngen))

        flags=cl.merger_flags(len(m1_ngen))

            
        BBH.m1=np.copy(np.array(m1_ngen))
        if(len(BBH.m1)==0):
            print("TDF too long: last gen is ", gen)
            print(BBH.m1)
            break
        BBH.idi=np.copy(np.array(idi_ngen))
        BBH.chi1=np.copy(np.array(chi1_ngen))
            
        tmerg=np.copy(np.array(tmerg_ngen))
        t3bb=np.copy(np.array(t3bb_ngen))
        
        Ngen=np.copy(np.array(Ngen))

        flags.flagSN=np.copy(np.array(flagSN_ngen))
        flags.flag_evap=np.copy(np.array(flag_evap_ngen))
        flags.flag_t3bb=np.copy(np.array(flag_t3bb_ngen))

        
        SC.rhocgs=np.copy(np.array(rhocgs_ngen))
        SC.sigma=np.copy(np.array(sigma_ngen))
        SC.maverage=np.copy(np.array(maverage_ngen))
        SC.vesc=np.copy(np.array(vesc_ngen))

        SC.Mtot=np.copy(np.array(Mtot_ngen))
        SC.SClifetime=np.copy(np.array(SClifetime_ngen))
        SC.feq=np.copy(np.array(feq_ngen))
        SC.binfrac=np.copy(np.array(binfrac_ngen))
        SC.rho=np.copy(np.array(rho_ngen))
        SC.trh0=np.copy(np.array(trh0_ngen))

        SC.cdensity=const.c2hm_dens * SC.rho/SC.maverage
        SC.sigma1D=SC.sigma/np.sqrt(3.)
        SC.rh=np.copy(np.array(rh_ngen))
        
        if(gf.evol_cluster==True):
            vesc0=np.copy(np.array(vesc0_ngen))
            rh0=np.copy(np.array(rh0_ngen))
            
            #SC.vesc,SC.sigma,SC.rh,SC.rho=auxi.update_vesc_rh_rho(vesc0,rh0,SC.trh0,SC.Mtot,tmerg+t3bb)
            SC.vesc,SC.sigma,SC.rh,SC.rho=auxi.nsc_prop(NC,(tmerg+t3bb))
            
            
            #print("sigma: ",SC.sigma)
            print("rho_ngen: ", SC.rho) 
            #print("vesc_ngen: ", SC.vesc/1e5)
            #print("time:",tmerg+t3bb)
            #print("rh_ngen:", SC.rh)
            SC.rhocgs=SC.rho * const.msun/(const.parsec**3.)
            #print(vesc_ngen)

        if gf.plotta == True:
            print(gen)
            plt.hist(t3bb)
            plt.xlabel('T$_\mathrm{DF}$ [Myr]')
            plt.show()
    

    
        tdel=np.zeros(len(BBH.m1))
        tmerg=np.copy(tmerg+t3bb)


        BBH = initngen.initialize_nthgen(BBH,SC,max1g,flags,gf) 


            
        if(gf.evol_cluster==True):
            BBH, SC, tdel, tmerg, flags = makegen.make_generation_evol(BBH, SC, tmerg, vesc0, rh0, mBHaverage, m1ii, flags, gf, NC)

        else:
            BBH, tdel, tmerg, flags = makegen.make_generation(BBH, SC, tmerg, mBHaverage, m1ii, flags, gf)


                    
        #print next generations
        if(gf.evol_cluster==True):
            for i in range(len(BBH.m1)):
                Ngen[i]=gen
                f2.write(str(BBH.idi[i])+" "+str(BBH.m1[i])+" "+str(BBH.m2[i])+" "+str(BBH.chi1[i])+" "+str(BBH.chi2[i])+" "+str(BBH.th1[i])+" "+str(BBH.th2[i])+" "+str(BBH.sma[i])+" "+str(BBH.ecc[i])+" "+str(t3bb[i])+" "+str(BBH.sma2[i])+" "+str(BBH.ecc2[i])+" "+str(tdel[i])+" "+str(tmerg[i])+" "+str(BBH.vk[i]/1e5)+" "+str(BBH.mmerg[i])+" "+str(BBH.chimerg[i])+" "+str(SC.vesc[i]/1e5)+" "+str(flags.flag[i])+" "+str(flags.flag2[i])+" "+str(flags.flag3[i])+" "+str(flags.flagSN[i])+" "+str(flags.flag_exch[i])+" "+str(flags.flag_t3bb[i])+" "+str(flags.flag_evap[i])+" "+str(SC.Mtot[i])+" "+str(BBH.ecc10[i])+" "+str(Ngen[i])+"\n")
                #if(flags.flag_t3bb[i]!="t3bb_ok"):
                #    print("horror t3bb too long\n")
        else:
            for i in range(len(BBH.m1)):
                Ngen[i]=gen
                f2.write(str(BBH.idi[i])+" "+str(BBH.m1[i])+" "+str(BBH.m2[i])+" "+str(BBH.chi1[i])+" "+str(BBH.chi2[i])+" "+str(BBH.th1[i])+" "+str(BBH.th2[i])+" "+str(BBH.sma[i])+" "+str(BBH.ecc[i])+" "+str(t3bb[i])+" "+str(BBH.sma2[i])+" "+str(BBH.ecc2[i])+" "+str(tdel[i])+" "+str(tmerg[i])+" "+str(BBH.vk[i]/1e5)+" "+str(BBH.mmerg[i])+" "+str(BBH.chimerg[i])+" "+str(SC.vesc[i]/1e5)+" "+str(flags.flag[i])+" "+str(flags.flag2[i])+" "+str(flags.flag3[i])+" "+str(flags.flagSN[i])+" "+str(flags.flag_exch[i])+" "+str(flags.flag_t3bb[i])+" "+str(flags.flag_evap[i])+" "+str(SC.Mtot[i])+" "+str(BBH.ecc10[i])+" "+str(Ngen[i])+"\n")
                #if(flags.flag_t3bb[i]!="t3bb_ok"):
                #    print("horror t3bb too long\n")                        

        if gf.plotta == True:
            plt.hist(BBH.vk/1e5)
            plt.xlabel("V$_\mathrm{k}$ [km/s]")
            plt.show()
            print(len(flags.flag2),len(flags.flagSN),len(flags.flag_t3bb))
            v=np.where((flags.flag2=="merg") & (flags.flag_t3bb=="t3bb_ok")) #all merged systems
            vv=v[00]
            print(flags.flag2[vv])

            plt.hist(BBH.mmerg[vv])
            plt.xlabel("M$_\mathrm{merg}$ [M$_\odot$]")
            plt.show()

            plt.hist(BBH.chimerg[vv])
            plt.xlabel("a$_\mathrm{merg}$")
            plt.show()
                
        
        gen+=1
    #closes else line 623 (if(len(a)=0)
#closes while(len(a)>0): line 515
ft.close()
f2.close()
