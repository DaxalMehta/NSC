import numpy as np
import auxiliary as auxi
import constants as const

import params as par
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

