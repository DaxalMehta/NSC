import numpy as np
import constants as const

###########################################
##############GENERAL PARAMETERS###########
###########################################

##############################################
### decide star cluster type (YSC, GC, NSC) ##
##############################################

star_cluster_type= "NSC"


##############################
### decide channel (Dyn, Orig)
#############################

channel = "Dyn"


###############################
### path of initial conditions
###############################

# common envelope parameter (A1 means alpha = 1)
alpha = "A1" 

# core collapse supernova model (rapid, delayed)
SN = "rapid" 

# directory of initial conditions
IC = "input/"+SN+"_"+alpha 


###################################
### number of 1g BBHs to simulate
###################################

BBH_number=1e4



##################################
### set progenitor's metallicity
##################################

Z = '0.02'


#########################################
### set rms of Maxwellian spin magnitude
### set a value in 0,1
#########################################

spin = 0.1


#########################################
### activate SN kick  vSN
### (if yes the BH is ejected when
### vSN>vesc)
#########################################

SN_kick = True


#########################################
### produce on the fly plots (True/False)
### option for debugging only
#########################################

plot = False

#########################################
### switch on cluster evolution (True/False)
### still WIP
#########################################

evol = True

#########################################
### switch on dynamical exchanges (True/False)
### still WIP
#########################################

exchanges = False



################################################
##############STAR CLUSTER PROPERTIES###########
################################################

userprovided = False
#### if userprovided = False FASTCLUSTER USES ITS OWN VALUES FOR SC PROPERTIES
#### if userprovided = True  FASTCLUSTER TAKES THE FOLLOWING VALUES

### SClifetime set up the star cluster lifetime
SClifetime = const.tHubble


#### binfrac is the fraction of surviving binaries in the cluster
binfrac = 1.0

#### Median_Mtot is the median value of the lognormal SC mass distribution
Median_Mtot=np.log10(2e5)

#### Sigma_Mtot is the standard deviation of the lognormal SC mass distribution
Sigma_Mtot=0.4


#### Median_rho is the median value of the lognormal distribution of SC density at the half-mass radius
Median_rho=np.log10(2e6)

#### Sigma_rho is the standard deviation of the lognormal distribution of SC density at the half-mass radius
Sigma_rho=0.4


#### average stellar mass in Msun
maverage = 1.0

#### Spitzer's parameter zeta. If zeta = 1 (<1) we (do not) have equipartition 
zeta = 0.1

################################################
##############BLACK HOLE PROPERTIES###########
################################################


userBH = False
#### if userBH = False FASTCLUSTER USES ITS OWN VALUES FOR BH PROPERTIES
#### if userBH = True  FASTCLUSTER TAKES THE FOLLOWING VALUES (user provided)

m2min = 3.0 #minimum mass of secondary BH
xi = 1.0    # hardening parameter  xi 
ki = 0.01   # eccentricity boost parameter kappa

################################################
##############GALAXY PROPERTIES###########
################################################

galaxy_type="Late" #Late or Early
galaxy_mass=1.E9
GC_mass=5.E5
slope=1.2
k=0.14
nsc_alpha=-0.67
beta=1.76
r=0.01		#in Kpc
den=1.2  #np.random.uniform(0,0.5) #Density profile slope
eccen=0.5
