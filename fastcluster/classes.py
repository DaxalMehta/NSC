import numpy as np

import auxiliary as auxi
import constants as const

import params as par

#import cla
#gl=cla.gal_properties()

#####################GENERAL FLAGS ######################################
class general_flags(object):
    'Flags to set up the simulation'

    def __init__(self):
        

        self.flagNSC = par.star_cluster_type #"NSC" "YSC" "GC"  #which star cluster flavor?

        self.channel= par.channel  #it can be Dyn or Orig
        
        self.alphaCE= par.alpha #value of alphaCE (fiducial is 1)
        
        self.SNmodel= par.SN #type of SN model (rapid or delayed)


        self.Z_list = par.Z #metallicity
        
        self.spinrms = par.spin #magnitude of spins

        self.Nref = par.BBH_number #number of BBHs

        self.in_dir= par.IC #MAIN INPUT DIR


        self.plotta= par.plot #does some test plot while running

        self.evol_cluster= par.evol #allow rh to expand and vesc, sigma to decrease

        self.SN_kick= par.SN_kick #marks BHs that are ejected and prevent second generation masses that should be ejected

        self.exchanges= par.exchanges #activates exchanges inside peters integration

        self.ngng=False #always draws secondary black hole from a nth-g black hole (if false second generation are nthg-nthg only if m2>mBHmax)



        if(abs(self.spinrms-0.1)<1e-3):
            self.flagspin = "_chi01"
        if(abs(self.spinrms-0.01)<1e-3):
            self.flagspin="_chi001"
        if(abs(self.spinrms-0.3)<1e-3):
            self.flagspin="_chi03"

        if self.exchanges:
            self.flagexch="_exch"
        else:
            self.flagexch=""
            
        if self.ngng:
            self.flagng="_2g2g"
        else:
            self.flagng=""
            
        if self.evol_cluster:
            self.flagclustev="_clusterevolv"
        else:
            self.flagclustev="_noclusterevolv"


        #self.out_top="output/" #main output directory
        self.out_top=self.SNmodel+"_"+self.alphaCE+"_"+self.flagNSC+self.flagspin+"_output"+self.flagclustev+self.flagexch+self.flagng+"/"

        self.run_name= self.channel+"/" 
        self.out_dir=self.out_top+self.run_name
        self.out_met=self.out_dir+self.Z_list
        

        
        #self.out_top=flagNSC+flagspin+"_output"+flagclustev+flagexch+flagng+"_sigma06"+"/" 


############STAR CLUSTER PROPERTIES######################
class SC_properties(object):
    'star cluster properties'
    def __init__(self, flagNSC, numBH):
        self.flagNSC = flagNSC
        self.numBH = numBH
        
        self.SClifetime=np.zeros(numBH,float)
        self.Mtot=np.zeros(numBH,float)

        self.rho=np.zeros(numBH,float)
        self.rhocgs=np.zeros(numBH,float)
        self.cdensity=np.zeros(numBH,float)
        
        self.vesc=np.zeros(numBH,float)

        self.sigma=np.zeros(numBH,float)
        self.sigma1D=np.zeros(numBH,float)
        
        self.maverage=np.zeros(numBH,float)

        
        self.feq=np.zeros(numBH,float)

        self.mmax=np.zeros(numBH,float)
    
        self.trh0=np.zeros(numBH,float)
        
        self.rh=np.zeros(numBH,float)
    
        self.binfrac=np.zeros(numBH,float) #hard binary fraction at equilibrium
    

    def init_SC(self, flagNSC,numBH,time): 
        for i in range(numBH):
            self.maverage[i] =1.0  #average star mass Msun
            self.feq[i] = 1.0 #equipartition fraction
            
            if(flagNSC=="NSC"):
                self.SClifetime[i] = const.tHubble #lifetime of the star cluster
                self.binfrac[i]=0.01
                #self.Mtot[i]=auxi.set_Mtot(np.log10(1.5e6),0.4)
                #self.rho[i]=auxi.set_rho(np.log10(1e5),0.4) 
                gl=gal_properties()
                self.Mtot[i]=gl.mass
                self.vesc[i],self.sigma[i],self.rh[i],self.rho[i]=auxi.nsc_prop(gl,time)
             
            elif(flagNSC=="GC"):
                self.SClifetime[i] = const.tHubble #was 12.0e3 #12 Gyr, lifetime of a Globular cluster (Gratton et al. 1997, 2003;VandenBerg et al. 2013) 
                self.binfrac[i]=0.1
                self.Mtot[i]=auxi.set_Mtot(np.log10(4e5),0.4) #was 2e5 then 5e5, then 4e5
                self.rho[i]=auxi.set_rho(np.log10(5e3),0.4) #was 2e3, then 1e4, then 4e3

            elif(flagNSC=="YSC"):
                self.SClifetime[i] = 1.0e3 #1 Gyr in Myr
                self.binfrac[i]=1.0
                self.Mtot[i]=auxi.set_Mtot(np.log10(2e4),0.4) 
                self.rho[i]=auxi.set_rho(np.log10(2e3),0.4) 


            if(par.userprovided==True): #read external properties
                self.SClifetime[i] = par.SClifetime
                self.binfrac[i]= par.binfrac
                self.Mtot[i]=auxi.set_Mtot(par.Median_Mtot,par.Sigma_Mtot) 
                self.rho[i]=auxi.set_rho(par.Median_rho,par.Sigma_rho) 
                self.maverage[i] = par.maverage  #average star mass Msun
                self.feq[i] = par.zeta #equipartition fraction
            

        #self.Mtot=10.**self.Mtot
        #self.rho=10.**self.rho
        #self.vesc = 40.0e5 * (self.Mtot/1e5)**(1./3.) * (self.rho/1e5)**(1./6.)

        self.cdensity= const.c2hm_dens * 1 * self.rho/self.maverage #central density pc^-3
        
        self.rhocgs=self.rho*const.msun/(const.parsec**3.)
        #Mtot = 1.0e5 * (vesc/40.0e5)**3. * np.sqrt(1.0e5/rho) #total star cluster mass Msun
        
        self.sigma1D = self.vesc/(2.*np.sqrt(3.)) #central 1D velocity dispersion  cm/s
        self.sigma = self.vesc/2. #central 3D velocity dispersion cm/s

        self.mmax=self.Mtot*1e-3 #maximum BH mass
        self.trh0=7.5*(self.Mtot/1e5)/(self.rho/1e5)**0.5 #calculates the initial half-mass relax time in Myr
        #self.rh=(3.*self.Mtot/(8.*np.pi*self.rho))**(1./3.)  #calculates the initial half mass radius in pc



        
    def get_SC(self,i):

        lista=[
            self.SClifetime[i],
            self.Mtot[i],

            self.rho[i],
            self.rhocgs[i],
            self.cdensity[i],
        
            self.vesc[i],

            self.sigma[i],
            self.sigma1D[i],
        
            self.maverage[i],

        
            self.feq[i],

            self.mmax[i],
    
            self.trh0[i],
        
            self.rh[i],
    
            self.binfrac[i]]


        return lista
    #return SClifetime,vesc,maverage,cdensity,feq,rho,rhocgs,Mtot,sigma1D,sigma,mmax,trh0,rh,binfrac



##############BBH properties#######################
class BBH_properties(object):
    'BBH properties'
    #m2=secondary (Msun)  sma=semimajor axis (Rsun) ecc=eccentricity chi1,chi2=dimensionless spin primary/secondary th1, th2=tilt angle (0 means aligned with orb.ang.mom)
    def __init__(self, numBH):
        self.numBH = numBH
        
        self.idi=np.zeros(numBH,int)
        
        self.m1=np.zeros(numBH,float)
        self.m2=np.zeros(numBH,float)
        
        self.sma=np.zeros(numBH,float)
        self.ecc=np.zeros(numBH,float)

        self.sma2=np.zeros(numBH,float)
        self.ecc2=np.zeros(numBH,float)

        self.ecc10=np.zeros(numBH,float) #eccentricity at fGW=10 Hz

        self.chi1=np.zeros(numBH,float)
        self.chi2=np.zeros(numBH,float)

        self.th1=np.zeros(numBH,float)
        self.th2=np.zeros(numBH,float)

        self.mmerg=np.zeros(numBH,float)
        self.chimerg=np.zeros(numBH,float)
        self.vk=np.zeros(numBH,float)
        
        self.m2min=3.0 #min BH mass in Msun
        self.xi=3.0 #constant for hardening 0.1-10 Hills 1983, Quinlan 1996, Sesana 2006
        self.ki=0.1 #constant for eccentricity variation 0.01-0.1 Quinlan 1996, Sesana 2006
        if par.userBH == True:
            self.m2min=par.m2min
            self.xi=par.xi
            self.ki=par.ki


#####################MERGER FLAGS ######################################
class merger_flags(object):
    'Flags to decide merger properties'

    def __init__(self, N):
        
        self.N = N
        self.flag_t3bb=np.empty(N,dtype=object)

        self.flagSN=np.empty(N,dtype=object)
            
        self.flag=np.empty(N,dtype=object)
        self.flag2=np.empty(N,dtype=object)
        self.flag3=np.empty(N,dtype=object)
        self.flag_evap=np.empty(N,dtype=object)

        self.hard=np.zeros(N,int) #zero means soft, one means hard
        
        self.flag_exch=np.zeros(N,int)
  
###################GALAXY PROPERTIES######################################
class gal_properties(object):
	def __init__(self):
		self.k=par.k
		self.slope=par.slope
		self.alpha=par.nsc_alpha
		self.beta=par.beta
		self.r=par.r
		self.den=par.den
		self.eccen=par.eccen
		self.gal_type=par.galaxy_type
		self.gmass=1e9
		self.gcmass=5.E5
		
		c1=2.**(1./(3.-self.slope))-1.
		self.Rg=2.37 * c1 * (self.gmass/1.E11)**self.k
		self.mass = (self.gmass * (self.r/(self.r + self.Rg))**(3.0-self.den))/(3.-self.den)
