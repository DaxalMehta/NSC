import numpy as np
import matplotlib.pyplot as plt
import classes 
from matplotlib import colors
from astropy import constants as const
from astropy import units as u
cl=classes.initial()

def radius(cl,i):
	mg=cl.gmass[i]/1.E11
	c1=2.**(1./(3. - cl.slope))-1.
	Rg=2.37 * c1 * (mg)**cl.k
	R=31.26 * Rg * mg**0.167
	return Rg,R

def density(cl,i):
	mg=cl.gmass[i]/1E11
	Rg,R=radius(cl,i)
	rho0=(3. - cl.den) *mg/(4. * np.pi * Rg**3.)
	rho=rho0/((cl.r/Rg)**cl.den * ((cl.r/Rg + 1.)**(4. - cl.den)))
	return rho0,rho	

def mass(cl,i):
        rho0,rho=density(cl,i)
        mass=rho * 4./3. * np.pi * cl.r**3
        ## the correct mass is
        Rg,R=radius(cl,i)
        mass = cl.gmass[i]/1.E11 * (cl.r/(cl.r + Rg))**(3.0 - cl.den)
        mass = mass/(3. - cl.den)
        return mass

def infall(cl,i):
        mg=cl.gmass[i]/1.E11
        mgc=cl.gcmass/1.E11
        a=((2.**(-3. + cl.slope))/(2.**(1./(3. - cl.slope)) - 1.0)**1.5)
        b=((mg)**(3./2. * (cl.k + 1.) - cl.alpha))
        c=((mgc)**(cl.alpha))

        
        a1 = 2.63
        a2 = 2.26
        a3 = 0.9        
        gfact = (2.-cl.slope)*((a1*pow(1./(2.-cl.slope), a2) + a3)*(1.-cl.eccen)+cl.eccen)

        f = 0.3 * gfact * pow(2.37,1.5)*(2.**(1./(3.-cl.slope))-1.0)**(1.5)

        rateg = 0.01*2.0**(-3.+cl.den)/f * mg**(1.5*(1.-cl.k)+cl.alpha) * mgc**(-1.-cl.alpha)

        #rateg=0.027 * a * b * c

        rateh=(2.0**(2.0-cl.slope))*(2.0**(1.0/(3.0-cl.slope))-1.0)**cl.beta*rateg
        return rateg,rateh


def extract(cl,i):
        mu = pow(np.random.random(), 1./(3.-cl.slope))
        Rg,R=radius(cl,i)
        pos = Rg * mu / (1.-mu)
        return pos

def timescales(cl,i):
        mg=cl.gmass[i]#/1e11
        mgc=cl.gcmass#/1e11
        rg,R=radius(cl,i)
        #rgc=np.random.random()*rg
        #the extraction must to be done according to the cluster density profile
        rgc=extract(cl, i)        

        rateg,rateh=infall(cl,i)
        time=rateg**-1.0

        #tburn=1.25e-2*2.0**cl.slope*mg*(mgc*rateg)**-1.0  ## a minus in the exponent is missing in the paper
        tburn = time * mg/mgc*0.01*pow(2.,-3.+cl.den)

        mg1=cl.gmass[i]/1e11
        mgc1=cl.gcmass/1e11

        #tdf=0.3*(2-cl.slope)*(4.93-3.93*cl.eccen)*(rg**3/mg1)**0.5*(mgc1/mg1)**cl.alpha*(rg/rg)**cl.beta
        #actual df time from Arca Sedda et al 2015
        a1 = 2.63
        a2 = 2.26
        a3 = 0.9
        rrh= rg / (pow(2.,1./(3.-cl.slope)) - 1.)
        gfact = (2-cl.slope)*((a1*pow(1./(2.-cl.slope), a2) + a3)*(1.-cl.eccen)+cl.eccen)
        tdf   = 0.3*(rg**3/mg1)**0.5* gfact * (mgc1/mg1)**cl.alpha * (rg/rg)**cl.beta

        print("warning: ",rg, mg1, mgc1, cl.beta, cl.alpha)
        return time,tburn,tdf


def escape_vel(cl,i,mass):
	rho=density(cl,i)[1]*1.E2
	#vesc=40.0E5 * (mass/1.E5)**(1./3.) * (mass/1.E8)**(1./6.)
	vesc=40.0E5 * (mass/1.E5)**(1./3.) * (rho/1.E5)**(1./6.)
	return vesc
	
N=len(cl.gmass)
M=len(cl.age)
MSC=np.zeros((N,M),float)
Vesc=np.zeros((N,M),float)
time=np.zeros(N,float)
tburn=np.zeros(N,float)
tdf=np.zeros(N,float)
for i in range(N):
	r1,r2=infall(cl,i)
	Rg,R=radius(cl,i)
	time[i],tburn[i],tdf[i]=timescales(cl,i)
	mg=cl.gmass[i]
	mgc=cl.gcmass
	ma=r2*mgc
	mi=mass(cl,i)*1e11
	MSC[i,:]=ma*cl.age+mi
	Vesc[i,:]=escape_vel(cl,i,MSC[i,:])
	print("Galaxy Mass: %.2e Msun" %cl.gmass[i])
	print("GC Mass: %.2e Msun" %cl.gcmass)
	print("Radius under which all GC mergers into NSC within th: %.2e Kpc" %Rg)
	print("Infall Rate h: %.2e 1/Myr" %r2)
	print("Infall Rate g: %.2e 1/Myr" %r1)
	print("Replenishment time : %.2e Myr" %time[i])
	print("Mass accretion rate: %.2e Msun/Myr" %ma)
	print("Mass inside 10pc initial and end: [%.2e %.2e] Msun" %(mi, MSC[i,M-1]))
	print("Central Density initial and end: [%.2e %.2e] Msun/pc3 " %(density(cl,i)[1]*100, MSC[i,M-1]/1000))
	print("Central Density initial and end: [%.2e %.2e] Msun/pc3 " %(mi/1000, MSC[i,M-1]/1000))
	print("Escape Velocity initial and end: [%.2e %.2e] cm/s" %(Vesc[i,0], Vesc[i,M-1]))
	print("#########################")
#	plt.plot(cl.age,MSC[i,:])
#	plt.xlabel('Time Myr')
#	plt.ylabel('Mass accreted Msun')
#	plt.plot(cl.age,Vesc[i,:])
#	plt.pause(1)
#"""
mg=cl.gmass
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
#"""

