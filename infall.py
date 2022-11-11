import numpy as np
import matplotlib.pyplot as plt
import classes 
from matplotlib import colors
from astropy import constants as const
from astropy import units as u
cl=classes.initial()

def radius(cl,i):
	mg=cl.gmass[i]/1e11
	c1=2**(1/(3-cl.slope))-1
	Rg=2.37*c1*(mg)**cl.k
	R=31.26*Rg*mg**0.167
	return Rg,R

def density(cl,i):
	mg=cl.gmass[i]/1e11
	Rg,R=radius(cl,i)
	rho0=(3-cl.den)*mg/(4*np.pi*Rg**3)
	rho=rho0/((cl.r/Rg)**cl.den*((cl.r/Rg+1)**(4-cl.den)))
	return rho0,rho	

def mass(cl,i):
	rho0,rho=density(cl,i)
	mass=rho*4/3*np.pi*cl.r**3
	return mass

def infall(cl,i):
	mg=cl.gmass[i]/1e11
	mgc=cl.gcmass/1e11
	a=((2**(-3+cl.slope))/(2**(1/(3-cl.slope))-1)**1.5)
	b=((mg)**(3/(2*(cl.k+1))-cl.alpha))
	c=((mgc)**(cl.alpha))
	rateg=0.027 * a * b * c
	rateh=(2**(2-cl.slope))*(2**(1/(3-cl.slope))-1)**cl.beta*rateg
	return rateg,rateh

def timescales(cl,i):
	mg=cl.gmass[i]#/1e11
	mgc=cl.gcmass#/1e11
	rg,R=radius(cl,i)
	rgc=np.random.random()*rg
	rateg,rateh=infall(cl,i)
	time=rateg**-1
	tburn=1.25e-2*2**cl.slope*mg*(mgc*rateg)**-1
	mg1=cl.gmass[i]/1e11
	mgc1=cl.gcmass/1e11
	tdf=0.3*(2-cl.slope)*(4.93-3.93*cl.eccen)*(rg**3/mg1)**0.5*(mgc1/mg1)**cl.alpha*(rg/rg)**cl.beta
	return time,tburn,tdf

N=len(cl.gmass)
M=len(cl.age)
MSC=np.zeros((N,M),float)
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

	print("Replenishment time: %.2e Myr" %time[i])
	print("Galaxy Mass: %.2e Msun" %cl.gmass[i])
	print("GC Mass: %.2e Msun" %cl.gcmass)
	print("Radius under which all GC mergers into NSC within th: %.2e Kpc" %Rg)
	print("Density inside 10pc: %.2e" %density(cl,i)[1])
	print("Mass inside 10pc: %.2e" %mi)
	print("Infall Rate g: %.2e 1/Myr" %r1)
	print("Replenishment time : %.2e Myr" %time[i])
	print("Infall Rate h: %.2e 1/Myr" %r2)
	print("Mass accretion rate: %.2e Msun/Myr" %ma)
	print("#########################")
	plt.plot(cl.age,MSC[i,:])
	plt.xlabel('Time Myr')
	plt.ylabel('Mass accreted Msun')
	plt.show()


mg=cl.gmass
plt.scatter(mg,time)
plt.yscale('log')
plt.ylabel('Time Myr')
plt.xscale('log')
plt.xlabel('Galaxy Mass')
plt.loglog(mg,tburn)
plt.loglog(mg,tdf)
plt.show()


