import numpy as np
import matplotlib.pyplot as plt
import os
import jimenez as ji
import constants as const

c5=const.c**5.
G3=const.G**3

def set_Mtot(mean,std): 
    #generates total star cluster mass from lognormal
    Mtot=np.random.normal(loc=mean,scale=std)
    return Mtot

def set_rho(mean,std): 
    #generates density at half-mass radius from lognormal
    rho=np.random.normal(loc=mean,scale=std)
    return rho
        
def update_vesc_rh_rho(vesc0,rh0,trh0,Mtot,t): 
    #evolves cluster properties
    #input vesc0 in cm/s, rh0 in pc, trh0 and t must be the same time unit, Mtot in Msun
    vesc=vesc0*(0.15*t/trh0+1.)**(-1./3.)
    rh=rh0*(0.15*t/trh0+1.)**(2./3.)
    rho=3.*Mtot/(8.*np.pi*rh**3)
    sigma=vesc/2.
    #print("Evolving rho: ", rho, rh, vesc, rh0, vesc0, Mtot, t) 
    return vesc,sigma,rh,rho

def make_exchange(m1,m2,a,e,m1ii): 
    #calculates exchange
    #m1,m2,a in cgs, m1ii in msun
    draw=np.random.randint(0,len(m1ii))
    m3=m1ii[draw]*const.msun
    flag=0
    if(m3>=m2):
        if(m3>m1):
            tt=m1
            m1=m3
            m2=tt
        else:
            m2=m3

        #the commented version is the experiment version (no much changes to justify it)
        #new_e=0.0
        #while(new_e>e):
        #    new_e=np.random.rand()
        #    new_e=np.sqrt(new_e)
        #e=new_e
        e=np.random.rand()
        e=np.sqrt(e)
        flag=1
    return m1,m2,a,e,flag


def t_12init(m1,maverage,sma,sigma,cdensity,binfrac):
    #exchange timescale from
    #miller & lauburg2009
    #input is msun for m1 and maverage 
    #cm for sma, cm/s for velocities
    #pc^-3 for central density

    sma=sma/const.AU

    t12=3.0e9 * (0.01/binfrac) * (1.0e6/cdensity) * (sigma/50.e5) * (12./(2.*maverage+m1)) * (1./sma)

    t12=t12/1e6 #from yr to Myr
    #print("sma hard = ", sma, "t12 = ", t12)

    return t12 #returns t12init (star star BH) in Myr


def t_12b(m1,m2,sma,sigma1D,rho,maverage,feq,binfrac,mBHaverage):
    #exchange timescale from antonini & rasio (2016)
    #after first generation, when BHs are segregated and Spitzer unstable
    #input is g for m1 and m2, msun for maverage and mBHaverage, 
    #cm for sma, cm/s for velocities
    #g/cm^3 for density
    mBH=mBHaverage #assume that the average BH mass is 10Msun
    ratio=mBH/maverage 

    m1=m1/const.msun
    m2=m2/const.msun
    sma=sma/const.AU
    density=rho/const.msun*(const.parsec)**3/maverage
    #print("Dens= ",density,"m1= ",m1)

    t12=3.0e8/feq * (0.01/binfrac) * (1.e6/density) * (sigma1D/30.0e5) * np.sqrt(10./ratio) * (10./(m1+m2)) * (1./sma)
    #print("T12= ", t12)
    t12=t12*const.yr #from yr to sec
    return t12 #returns t12 in sec


def t_dis(m1,m2,sma,maverage,sigma,density):
    #disruption time of a binary (binney & tremaine 1987) no longer used
    density=density/(const.parsec**3)
    maverage=maverage*const.msun
    m1=m1*const.msun
    m2=m2*const.msun
    sma=sma*const.Rsun
    Eb=const.G*m1*m2/2./sma
    tdis=9.*Eb*Eb/(16. * np.sqrt(np.pi)*density*const.G*const.G* maverage**4 * sigma) * (1.+4.*maverage *sigma *sigma/(15.*Eb)) * (1.+np.exp(3.*Eb/(4.*maverage*sigma*sigma)))
    return tdis

    
def t_dynfric(sigma,m,rhocgs):
    #dynamical friction timescale from Chandrasekhar 1943
    #input sigma in cm/s, m in msun, rhocgs in g/cm^3
    coulomb_log=10.0
    tDF=3./(4.*(2.*np.pi)**0.5 * const.G * const.G * coulomb_log) * sigma**3/(m * const.msun * rhocgs)
    return tDF


def t_dynfric2(m1,maverage,Mtot,rh):
    #dynamical friction timescale from Portegies Zwart 2010
    #input sigma in cm/s, m in msun, rhocgs in g/cm^3
    tDF=200.*np.sqrt(Mtot/1e6) * (1.7/1.3 * rh/1.0)**1.5 * (maverage/m1)
    return tDF

def t_3bb(density,feq,sigma1D,m1):
    #time for formation of a binary by 3body encounters (3 single bodies)
    #lee 1995
    #input density in pc^-3, sigma1D in cm/s and m1 in Msun
    #fragione t3bb=125.0*(1e6/density)**2  * (sigma1D/30.e5/feq)**9 * (20./m1)**5 #time for three-body capture in Myr, Lee (1995)
    #antonini
    t3bb=125.0*(1e6/density)**2  * (sigma1D/30.e5/feq)**9 * (20./m1)**5 #time for three-body capture in Myr, Lee (1995)
    #t3bb=4.0e3*(1e6/density)**2  * (sigma1D/30.e5/feq)**9 * (20./m1)**5 #time for three-body capture in Myr, Lee (1995)
    return t3bb #returns t3bb in Myr

def a_ej(xi,maverage,m1,m2,vesc):
    #semi-major axis below which binary is ejected by 3body encounters
    #input masses in msun, vesc in cm/s
    maverage=maverage*const.msun
    m1=m1*const.msun
    m2=m2*const.msun
    aej=xi*maverage*maverage/(m1+m2)**3. * const.G * m1 * m2/vesc/vesc #semi-major axis for dynamical ejection MM: fixed bug was 2*xi / 2021/05/31 
    #2.*xi*maverage*maverage/(m1+m2)**3. * const.G * m1 * m2/vesc/vesc #semi-major axis for dynamical ejection 
    return aej

def a_GW(ecc,m1,m2,xi,sigma,rho):
    #semi-major axis below which GWs dominate evolution
    #obtained assuming adot_ej=adot_gw (Baibhav et al. 2020)
    #input masses in msun
    m1=m1*const.msun
    m2=m2*const.msun
    #rho= #must be central density in cgs
    #tH=const.tHubble*const.yr*1.0e6
    #aGW=(tH* 256./5.*G3/c5/(1.-ecc*ecc)**3.5 * m1 * m2 * (m1+m2))**(1./4.)
    aGW=32.* const.G**2./(5. * np.pi * xi * const.c**5.) * (sigma/rho) * m1 * m2 * (m1+m2)/((1-ecc*ecc)**3.5) * (1.+73./24.*ecc*ecc+37./96. * ecc * ecc * ecc * ecc)
    aGW=aGW**(1./5.)
    return aGW

#def semimajor_hardonly(m1,m2,maverage,sigma): 
#    #generates 1g semi-major axis
#    #now P(a)=1/a arbitrary
#    #input masses in msun,sigma in cm/s
#    lower=1e0 #rsun
#    upper=1e3 #rsun
#    maxi=np.log(upper)
#    mini=np.log(lower)
#    x=np.random.rand()
#    sma=np.exp(x*(maxi-mini)+mini)
#    flag="hard"
#    flag=check_hard(m1,m2,sma,maverage,sigma)
#    while(flag=="soft"):
#        x=np.random.rand()
#        sma=np.exp(x*(maxi-mini)+mini)
#        flag=check_hard(m1,m2,sma,maverage,sigma)
#    return sma


#def check_hard_hardonly(m1,m2,sma,maverage,sigma):
#    #check if a binary is hard (survives) or soft (breaks)
#    #input masses in msun, sma in Rsun, sigma in cm/s
#    Eb=const.G*m1*m2*const.msun*const.msun/(2.*sma*const.Rsun)
#    Ek=0.5*maverage*const.msun*sigma*sigma
#    flag="hard"
#    if(Eb<Ek):
#        #print("Warning: binary is soft")
#        flag="soft"    
#    return flag


def semimajor(m1,m2,lower=1e0,upper=1e3): 
    #generates 1g semi-major axis
    #now P(a)=1/a arbitrary
    #input masses in msun,sigma in cm/s
    #lower=1e0 #rsun
    #upper=1e3 #rsun
    N=len(m1)
    maxi=np.log(upper)
    mini=np.log(lower)
    x=np.random.rand(N)
    sma=np.exp(x*(maxi-mini)+mini)
    return sma


def check_hard(m1,m2,sma,maverage,sigma):
    #check if a binary is hard (survives) or soft (breaks)
    #input masses in msun, sma in Rsun, sigma in cm/s
    Eb=const.G*m1*m2*const.msun*const.msun/(2.*sma*const.Rsun)
    Ek=0.5*maverage*const.msun*sigma*sigma

    flag=np.zeros(len(Eb),dtype=int) #zeros for soft, ones for hard
    a=np.where(Eb>=Ek)
    #print("Warning: binary is hard")
    flag[a[0]]=np.ones(len(a[0]),dtype=int)
    print(len(flag)-len(flag[a[0]]))
    return flag


def ecc(m1): 
    #generates initial eccentricity
    #now P(e)=2 e arbitrary
    N=len(m1)
    x=np.random.rand(N)
    e=np.sqrt(x)
    return e


###################generate secondary (o'leary 2016)###############
def secondary(m1,m2min):
    #generates secondary mass (o'leary 2016)
    #N=len(m1)
    C=(2.*m1)**5.-(m1+m2min)**5.
    if(C<1e-5):
        print("horror C is zero", C)
    C=1./C
    y=np.random.rand()
    m2=(y/C+(m1+m2min)**5.)**(1./5.)-m1
    return m2


#################spin of the merger remnant ###############################
#generates final spin of the merger remnant from jimenez-forteza et al 2017
def final_spin(chi1,chi2,m1,m2):
    #input masses in msun
    if(m1<m2):
        print("m1 <m2 I FLIP ", m1, m2)
        m1,m2=m2,m1
        chi1,chi2=chi2,chi1

    chif = ji.bbh_final_spin_non_precessing_UIB2016(m1,m2,chi1,chi2,"v2")
    if(chif>1.0):
        print("horror= ", chif)

    return chif



#################spin of the merger remnant, rezzolla #######################
#final spin from rezzolla et al 2008 PhRv 2008PhRvD..78d4002R only for aligned system
#should be generalized using Hofmann et al. 2016
def final_spin_rezz(a1,a2,m1,m2):
    #input masses in msun
    if(m1<m2):
        print("m1 <m2 I FLIP ", m1, m2)
        m1,m2=m2,m1
        chi1,chi2=chi2,chi1

    nu=m1 * m2 /(m1+m2)/(m1+m2)
    q=m2/m1
    a=(a1+a2*q*q)/(1+q*q)

    s4=-0.129
    s5=-0.384
    t0=-2.686
    t2=-3.454
    t3=2.353
    afin=a + s4 * a* a * nu + s5 * a * nu * nu + t0 * a * nu + 2. * (3)**0.5 * nu +t2 * nu * nu + t3 * nu*nu*nu
    #if(afin>1.0):
    #    print("spin_final = ", afin)
    afin=min(0.998,afin)

    return afin


################# mass of the merger remnant###############################
#generates final mass of the merger remnant from jimenez-forteza et al 2017
def final_mass(m1,m2,chi1,chi2):
    #input masses in msun
    if(m1<m2):
        print("m1 <m2 I FLIP ", m1, m2)
        m1,m2=m2,m1
        chi1,chi2=chi2,chi1
    mf = ji.bbh_final_mass_non_precessing_UIB2016(m1,m2,chi1,chi2,"v2")
    if(mf>(m1+m2)):
        print("horror mf= ", mf, m1+m2)
    return mf




########################## SPIN MAGNITUDE ################################
# generates spin magnitude from a maxwellian
def spin(rms=0.1):
    a=2.0
    #print(rms)
    while(a>1.0):
        #default_scale =0.1
        #alternative_scale=0.01
        n1=np.random.normal(loc=0.0,scale=rms)
        n2=np.random.normal(loc=0.0,scale=rms)
        n3=np.random.normal(loc=0.0,scale=rms)

        a=(n1*n1+n2*n2+n3*n3)**0.5
    return a

########################## SPIN TILT ################################
def set_theta(x): #generate spin tilt theta
    theta=x*2.-1.0
    theta=np.arccos(theta)
    return theta
    
######################### RELATIVISTIC KICKS ################################
def relativistic_kick(m1,m2,a1,a2,th1,th2,version="maggiore"):

    #input masses in msun
    #parallel is perpendicular to orbital plane
    #perpendicular 1 and 2 are in the orbital plane
    a1per=a1*np.sin(th1)
    a1par=a1*np.cos(th1)

    a2per=a2*np.sin(th2)
    a2par=a2*np.cos(th2)
    
    phi1=np.random.rand()*2.*np.pi



    if((version=="campanelli") or (version=="lousto1") or  (version=="lousto2")):
        if(m1>m2):
            #print("m1 <m2 I FLIP ", m1, m2)
            m1,m2=m2,m1 #lousto use m1>m2
            a1,a2=a2,a1 #lousto use m1>m2

        q=m1/m2 #q=m1/m2 with m1<m2 Lousto & Zlochower 2011

    if(version=="maggiore"):
        if(m1<m2):
            print("m1 <m2 I FLIP ", m1, m2)
            m1,m2=m2,m1 #lousto use m1>m2
            a1,a2=a2,a1 #lousto use m1>m2

        q=m2/m1 #astrophysicists use m1>m2


    eta=q/(1.+q)**2



    Am=1.2e4 #A
    Bm=-0.93  #B
    H=6.9e3
    K=6.0e4

    
    x=145./360.*2.*np.pi




    ###VERSION CAMPANELLI ET AL 2007
    if(version=="campanelli"):

        vm = Am * eta * eta * (1.-4.*eta)**0.5 * (1.+Bm*eta)
        vper= H * eta * eta/(1.+q) * (a2par-q*a1par)
        vpar= K * eta * eta/(1.+q) * (a2per-q*a1per) * np.cos(phi1)

    ###VERSION LOUSTO ET AL. 2012 eq 4
    if(version=="lousto1"):

        HS=0.0 # note lousto never writes it down explicitly
        KS=-4.254
        BK=0.0
        BH=0.0


        phi2=np.random.rand()*2.*np.pi

        vm= Am * eta * eta * (1.-q)/(1.+q) * (1. + Bm * eta)

        vper = H * eta * eta / (1.+q) * ((1. + BH * eta) * (a2par - q * a1par) + HS * (1-q)/(1+q)**2 * (a2par + q * q * a1par))

        vpar = K * eta *eta /(1.+q) * ((1. + BK * eta) * abs(a2per - q * a1per) * np.cos(phi1) + KS * (1. - q)/(1. + q)**2 * abs(a2per + q * q * a1per) * np.cos(phi2))

    
    ###VERSION LOUSTO ET AL. 2012 eq 12
    if(version=="lousto2"):

        V11=3678.0 #lousto et al. 2012
        VA=2.481e3
        VB=1.792e3
        VC=1.507e3
        xx=2. * (a2par + q * q * a1par) / (1.+q)**2


        vm= Am * eta * eta * (1.-q)/(1.+q) * (1. + Bm * eta)

        vper = H * eta * eta / (1.+q) * (a2par - q * a1par) 

        vpar= 16. * eta * eta /(1.+q) * (V11 + VA * xx + VB * xx * xx + VC * xx * xx * xx) * abs(a2per - q * a1per)*np.cos(phi1)
    

     ###VERSION MAGGIORE BOOK GW, 2 eqs 14.203,4,5
     #note eta := q/(1+q)**2 --> q**2/(1+q)**5 = eta * eta /(1+q)**5

    if(version=="maggiore"):

        Cmagg = 457.0 #km/s
        Dmagg = 3750.0 #km/s

        vm= Am  * eta * eta * (1.-q)/(1.+q) * (1. + Bm * eta)

        vper = Cmagg * 16. * eta * eta / (1. + q) * abs(a1par - q * a2par) 

        vpar = Dmagg * 16. * eta * eta / (1. + q) * abs(a1per - q * a2per) * np.cos(phi1)

    vk=(vm*vm+vper*vper+2.*vm*vper*np.cos(x)+vpar*vpar)**0.5

    #print(m1,m2,a1,a2,th1,th2,vk)

    if((version=="campanelli") or (version=="lousto1") or  (version=="lousto2")):
        th1,th2=th2,th1

    return(1e5*vk)




######################## SN kick ##################################
#generates SN kick based on linear momentum conservation
def SN_kick(mbh):
    #input masses in msun
    mNS=1.33 #typical NS mass in Msun
    vg=np.zeros(3,float)
    for i in range(len(vg)):
        vg[i]=np.random.normal(loc=0.0,scale=265.) 
    v=np.linalg.norm(vg) #kick magnitude of a NS from Hobbs et al. 2005
    v=v*1e5 *mNS/mbh #kick modulated by linear momentum in cm/s
    return v
    
###################### GC infall ###################################
def infall(cl):
	mg=cl.gmass/1.E11
	mgc=cl.gcmass/1.E11
	a1 = 2.63
	a2 = 2.26
	a3 = 0.9        
	gfact=(2.-cl.slope) * ((a1 * pow(1./(2.-cl.slope), a2) + a3) * (1.-cl.eccen) + cl.eccen)
	f=0.3 * gfact * pow(2.37,1.5)*(2.**(1./(3.-cl.slope))-1.0)**(1.5)
	rateg=0.01*2.0**(-3. + cl.den)/f * mg**(1.5*(1.-cl.k) + cl.alpha) * mgc**(-1.-cl.alpha)
	rateh=(2.0**(2.0-cl.slope)) * (2.0**(1.0/(3.0-cl.slope))-1.0)**cl.beta*rateg
	return rateg,rateh


def extract(cl):
	mu = pow(np.random.random(), 1./(3.-cl.slope))
#	Rg,R=radius(cl)
	Rg=cl.Rg
	pos = Rg * mu/(1.-mu)
	return pos

def timescales(cl):
#	rg,R=radius(cl)
	rg=cl.Rg
	rgc=extract(cl)        #the extraction must to be done according to the cluster density profile
	rateg,rateh=infall(cl)
	time=rateg**-1.0
	mg=cl.gmass
	mgc=cl.gcmass
	#tburn=1.25e-2*2.0**cl.slope*mg*(mgc*rateg)**-1.0  ## a minus in the exponent is missing in the paper
	tburn =10 * time * mg/mgc*0.01*pow(2.,-3.+cl.den)

	mg/=1.E11
	mgc/=1.E11
	#tdf=0.3*(2-cl.slope)*(4.93-3.93*cl.eccen)*(rg**3/mg1)**0.5*(mgc1/mg1)**cl.alpha*(rg/rg)**cl.beta
	#actual df time from Arca Sedda et al 2015
	a1 = 2.63
	a2 = 2.26
	a3 = 0.9
	rrh= rg/(pow(2.,1./(3.-cl.slope))-1.)
	gfact = (2-cl.slope)*((a1*pow(1./(2.-cl.slope), a2) + a3)*(1.-cl.eccen) + cl.eccen)
	tdf   = 0.3*(rg**3/mg)**0.5* gfact * (mgc/mg)**cl.alpha * (rg/rg)**cl.beta

#	print("warning: ",rg, mg1, mgc1, cl.beta, cl.alpha)
	return time,tburn,tdf


def halfmass_rad(cl,mass):
	if cl.gal_type=="Early":
		c1,c2=6.27,1.95E6
		alpha,beta=0.347,-0.024
		reff=c1 * 10**(alpha * np.log10(mass/c2) + beta)
	elif cl.gal_type=="Late":
		c1,c2=3.31,3.65E6
		alpha,beta=0.321,-0.011
		reff=c1 * 10.**(alpha * np.log10(mass/c2) + beta)
	rhm=4/3*reff
	return rhm,reff

def scale_den(cl,mass):
	rhm=halfmass_rad(cl,mass)[1]
	den=(3. * mass)/(4. * np.pi * rhm**3)
	#den=mass/rhm**3.
	return den

def escape_vel(cl,mass):
	rho=scale_den(cl,mass)
	vesc=40.0E5 * (mass/1.E5)**(1./3.) * (rho/1.E5)**(1./6.)
	return vesc


def nsc_prop(cl,time):
	t=min(const.tHubble,timescales(cl)[1])
	time=np.where(time<t,time,t)
#	if time>timescales(cl)[1]: time=timescales(cl)[1]
	init_mass=cl.mass  			# MSun
	infall_rate=infall(cl)[1] 		# 1/Myr
	acc_rate=infall_rate*cl.gcmass 		#Msun/Myr
	NSC_mass=init_mass+acc_rate*time 	#Msun	NSC_mass=np.where(NSC_mass<1.E8,NSC_mass,1.E8)
	rhm=halfmass_rad(cl,NSC_mass)[1] 	#pc
	scale_d=scale_den(cl,NSC_mass) 		#Msun/pc3
	esc_vel=escape_vel(cl,NSC_mass) 	#cm/s
	v_disp=esc_vel/2.
	return esc_vel,v_disp,rhm,scale_d
  
