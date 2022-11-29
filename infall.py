import numpy as np
import matplotlib.pyplot as plt
import classes 
from matplotlib import colors
from astropy import constants as const
from astropy import units as u
cl=classes.initial()

def infall(cl):
	##Calculate the infall rate of globular clusters
	## Takes GC mass in msun, Galaxy mass in msun, several internal constants of the galaxy.
	## Output is in 1/Myr
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
	### Calculate the radius of GC according to its distribution
	## Takes Rg in pc, the radius under which all GC merge into NSC within hubble time.
	## Output in pc
	mu = pow(np.random.random(), 1./(3.-cl.slope))
#	Rg,R=radius(cl)
	Rg=cl.Rg
	pos = Rg * mu/(1.-mu)
	return pos

def timescales(cl):
	## Calculate the timescales it will take for th GCs to merge
	## Takes Rg in pc, infall rate in 1/Myr, GC mass and Galaxy mass in msun
	## Returns timescales in Myr
	rg=cl.Rg
	rgc=extract(cl)        #the extraction must to be done according to the cluster density profile
	rateg,rateh=infall(cl)
	time=rateg**-1.0
	mg=cl.gmass
	mgc=cl.gcmass
	#tburn=1.25e-2*2.0**cl.slope*mg*(mgc*rateg)**-1.0  ## a minus in the exponent is missing in the paper
	tburn = time * mg/mgc*0.01*pow(2.,-3.+cl.den)#

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
print(timescales(cl)[1])
def halfmass_rad(cl,mass):
	## Scaling relations between reff and Mnsc
	## Takes NSC mass as input in msun
	## Outputs reff in pc
	if cl.gal_type=="Early":
		c1,c2=6.27,1.95E6
		alpha,beta=0.347,-0.024
		reff=c1 * 10**(alpha * np.log10(mass/c2) + beta)
	elif cl.gal_type=="Late":
		c1,c2=3.31,3.60E6
		alpha,beta=0.321,-0.011
		reff=c1 * 10.**(alpha * np.log10(mass/c2) + beta)
	rhm=4/3*reff	
	return rhm,reff


def scale_den(cl,mass):
	## Density of NSC at half-mass radius
	## Takes NSC mass in msun and reff in pc
	## Returns density in msun/pc3
	
	rhm=halfmass_rad(cl,mass)[1]
	#den=(3. * mass)/(4. * np.pi * rhm**3)
	den=(mass)/(rhm**3)
	return den
def escape_vel(cl,mass):
	## Escape velocity from NSC
	## Takes density in msun.pc3
	## Returns velocity in cm/s
	
	rho=scale_den(cl,mass)
	vesc=40.0E5 * (mass/1.E5)**(1./3.) * (rho/1.E5)**(1./6.)
	return vesc


def nsc_prop(cl,time):
	init_mass=cl.mass
	infall_rate=infall(cl)[1] 		# 1/Myr
	acc_rate=infall_rate*cl.gcmass 		#Msun/Myr
	NSC_mass=init_mass+acc_rate*time 	#Msun
	
	rhm=halfmass_rad(cl,NSC_mass)[1] 	#pc	
	scale_d=scale_den(cl,NSC_mass) 		#Msun/pc3
	esc_vel=escape_vel(cl,NSC_mass) 	#cm/s
	
	return NSC_mass,rhm,scale_d,esc_vel

print("%.3e %.3e %.3e %.3e" %(nsc_prop(cl,1000)))


"""


mg=cl.gmass
plt.loglog(mg,MSC[:,M-1])
plt.xlabel('Log(Mg)')
plt.ylabel('Log(Msc)')
plt.xlim(1.E6,1.E11)

plt.show()
plt.scatter(mg,time)
plt.yscale('log')
plt.ylabel('Time Myr')
plt.xscale('log')
plt.xlabel('Galaxy Mass')
plt.ylim(0.1,1.E5)
plt.grid(which="both")
plt.loglog(mg,tburn,c="green")
plt.loglog(mg,tdf,c="red")
plt.savefig("out.jpeg")
plt.close()
"""

