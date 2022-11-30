FASTCLUSTER by Michela Mapelli: calculates binary black hole mergers in young, globular and nuclear star clusters following hierarchical mergers
-------------------------------------------------------------------------------

FASTCLUSTER is an open-source population-synthesis code for binary black hole dynamics in star clusters. You can use FASTCLUSTER to study the hierarchical formation of binary black holes. We hope you will enjoy it. 


RELATED PAPERS:
--------------
We used FASTCLUSTER to produce the results of the following papers:


 Mapelli et al. 2021, MNRAS, 505, 339, https://ui.adsabs.harvard.edu/abs/2021MNRAS.505..339M/abstract

Mapelli et al. 2022, MNRAS, in press, https://arxiv.org/abs/2109.06222


You can have a look at them for details about the code. Please, do not forget to add a citation to our papers if you use FASTCLUSTER for a publication.

For questions: michela.mapelli@unipd.it


RUN
----

python fastcluster.py


FILES
------

**fastcluster.py:**  main, calculates evolution of dynamical and original binaries

**params.py:**	 parameter file, contains the user-provided input parameters


auxiliary.py:	  contains auxiliary functions

classes.py:	 defines fundamental classes:
		 general_flags 
		 SC_properties
		 BBH_properties
		 merger_flags

constants.py:	  contains main constants

initialize_dynbin.py:	initializes first generation dynamical binary

initialize_origbin.py:	initialize original binary

initialize_nthgen.py:	initialize nth-generation binary


jimenez.py:	  code from Jimenez-Forteza et al. (2017, DOI: 10.1103/PhysRevD.95.064024)
		  to calculate mass and spin of remnant,
		  arxiv: https://arxiv.org/abs/1611.00332

make_generation.py:	main calculation with Peters and dynamical encounters
			(calls peters/peters_exchanges),
			set up the relevant flags,
			if merger calculates remnant mass/spin
			and relativistic kick velocity


peters.py:	  main integrator of semi-major axis and eccentricity evolution


analysis/	  contains scripts for the analysis of the output



OUTPUTS
-------

**1) file "first_gen.txt": each line contains the information for one of the first-generation BBHs**
---------------------------------------------------------

col 0:identifier of the BBH (a number from 1 to the maximum number of BBHs)

col 1:M1/Msun, mass of the primary black hole in Msun

col 2:M2/Msun, mass of the secondary black hole in Msun

col 3:chi1, dimensionless spin of the primary black hole

col 4:chi2, dimensionless spin of the secondary black hole

col 5:theta1, tilt angle between the orbital angular momentum of the binary and the spin of the primary black hole

col 6:theta2, tilt angle between the orbital angular momentum of the binary and the spin of the secondary black hole

col 7:SMA(Rsun), initial semi-major axis of the BBH in solar radii

col 8:ecc, initial eccentricity of the binary system

col 9:tDF+min(t12,t3bb)/Myr, 
    formation time of the first generation BBH in Myr


col 10:SMAfin(cm), final semimajor axis in cm

col 11:eccfin, final eccentricity

col 12:tpeters/Myr, time of gravitational wave decay+hardening of the binary system in Myr

col 13:(tDF+t3bb+tpeters)/Myr, total time from the formation to the merger of the BBH in Myr (delay time)

col 14:vkick/kms, relativistic kick in km/s

col 15: mrem/Msun, mass of the compact remnant in solar masses

col 16:arem, spin of the compact remnant

col 17:vesc/kms, escape velocity from the stellar cluster (in km/s)

col 18:flag1, internal flag

col 19:flag2, internal flag

col 20:flag3, internal flag

col 21:flagSN, internal flag

col 22:flag_exch, internal flag

col 23:flag_t3bb, internal flag

col 24:flag_evap, internal flag

col 25:Mtot/Msun, total mass of the star cluster in solar masses

col 26:ecc(10Hz), eccentricity of the binary system at 10 Hz

col 27:Ngen, generation number




**2) file "nth_generation.txt": same as first_gen.txt, but for Nth generation BBHs. Each line is a single BBH.**
-----------------------------------------------------------------------------

col 0:identifier of the BBH (the identifier is the same as the first generation BBH in each generation, to track them from one generation to the other)

col 1:M1/Msun, mass of the primary black hole in Msun

col 2:M2/Msun, mass of the secondary black hole in Msun

col 3:chi1, dimensionless spin of the primary black hole

col 4:chi2, dimensionless spin of the secondary black hole

col 5:theta1, tilt angle between the orbital angular momentum of the binary and the spin of the primary black hole

col 6:theta2, tilt angle between the orbital angular momentum of the binary and the spin of the secondary black hole

col 7:SMA(Rsun), initial semi-major axis of the Ng BBH in solar radii

col 8:ecc, initial eccentricity of the Ng binary system

col 9:tDF+min(t12,t3bb)/Myr, 
    formation time of the Ng generation BBH in Myr,
    same as for 1st generation

col 10:SMAfin(cm), final semimajor axis in cm

col 11:eccfin, final eccentricity

col 12:tpeters/Myr, time of gravitational wave decay+hardening of the binary system in Myr

col 13:(ngen (tDF+t3bb+tpeters))/Myr, total time from the formation of the 1st generation BBH to the merger of the Ng BBH in Myr (delay time)

col 14:vkick/kms, relativistic kick in km/s

col 15: mrem/Msun, mass of the compact remnant in solar masses

col 16:arem, spin of the compact remnant

col 17:vesc/kms, escape velocity from the stellar cluster (in km/s)

col 18:flag1, internal flag

col 19:flag2, internal flag

col 20:flag3, internal flag

col 21:flagSN, internal flag

col 22:flag_exch, internal flag

col 23:flag_t3bb, internal flag

col 24:flag_evap, internal flag

col 25:Mtot/Msun, total mass of the star cluster in solar masses

col 26:ecc(10Hz), eccentricity of the binary system at 10 Hz

col 27:Ngen, generation number



**3) timescales_1stgen.txt**
-------------------------
some timescales of 1st generation BBHs

col 0: identifier of BBH

col 1: t_DF, dynamical friction timescale (Myr)

col 2: t_3bb, 3body capture timescale (Myr)

col 3: t_12, exchange timescale (Myr)

col 4: t_form, time of supernova explosion (Myr)


**4) timescales_nthgen.txt**
--------------------------
same as timescales_1stgen.txt but for Ng BBHs

col 0: identifier of BBH

col 1: t_DF, dynamical friction timescale (Myr)

col 2: t_3bb, 3body capture timescale (Myr)

col 3: t_12, exchange timescale (Myr)

