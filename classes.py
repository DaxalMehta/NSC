import numpy as np
import inputs
class initial(object):
	def __init__(self):
		self.gmass=inputs.gmass
		self.k=inputs.k
		self.slope=inputs.slope
		self.alpha=inputs.alpha
		self.beta=inputs.beta
		self.r=inputs.r
		self.den=inputs.den
		self.eccen=inputs.eccen
		self.gal_type=inputs.galaxy_type
		self.gcmass=inputs.GC_mass
		c1=2.**(1./(3.-self.slope))-1.
		self.Rg=2.37 * c1 * (self.gmass/1.E11)**self.k #Radius under which all GC fall into NSC in hubble time
		self.mass = (self.gmass * (self.r/(self.r + self.Rg))**(3.0-self.den))/(3.-self.den)
		#rho0=(3.-self.den) * (self.gmass/1.E11)/(4. * np.pi * self.Rg**3.)
		#self.rho=rho0/((self.r/self.Rg)**self.den * ((self.r/self.Rg + 1.)**(4.-self.den)))
