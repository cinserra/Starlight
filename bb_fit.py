#! /usr/bin/python

# Created by Cosimo Inserra
# GitHub: https://github.com/cinserra
# @COSMO_83
from scipy import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate

################### for the help ##################
from optparse import OptionParser

description = "Script to evaluate the blackbody temperature and spectral energy distribution of astronomical sources"
usage = "%prog "
if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description,)
    parser.add_option("-v", "--verbose",dest="verbose",\
                  action="store_true",default=False,
                  help='Print tasks description')
    parser.add_option("-m", "--magback",dest="magback", action="store",default='', type="str",
                  help='Do you want to retrieve some magnitudes? [yes/no; default=no]')
    option,args = parser.parse_args()


########### Options that can be changed thanks to option parser #########
if option.magback == 'yes':
    _magback = 'yes'
else:
    _magback = 'no'

'''
# example usage! #

bands = 'griz'

wl = []
fref = []

for i in bands:
    wl.append(wle[i])
    fref.append(zp[i]*1e-11) # because reference fluxes in units of 10^-11 erg/s/cm2/A

wl = np.array(wl)
fref = np.array(fref)

# note: can struggle to fit u/U sometimes, since SN isn't a good BB in the UV!

mags = np.array([19.34, 19.65, 19.88, 19.89])

z = 0.26

# cosmocalc - for getting SN distance from redshift
# Credits: Ned Wtight and James Schombert for the python adaptation
# These are the constant the user should change (if a different cosmology is used)

H0 = 72                         # Hubble constant
WM = 0.27                        # Omega(matter)
WV = 1.0 - WM - 0.4165/(H0*H0)  # Omega(vacuum) or lambda

WR = 0.        # Omega(radiation)
WK = 0.        # Omega curvaturve = 1-Omega(total)
c = 299792.458 # velocity of light in km/sec
Tyr = 977.8    # coefficent for converting 1/H into Gyr
DTT = 0.0      # time from z to now in units of 1/H0
DCMR = 0.0     # comoving radial distance in units of c/H0
DA = 0.0       # angular size distance
DL = 0.0       # luminosity distance
DL_Mpc = 0.0
a = 1.0        # 1/(1+z), the scale factor of the Universe
az = 0.5       # 1/(1+z(object))
 .....
 .....

# end of example parameters #
'''

### definitions and lists ####
def bbody(lam,T,A):
    Blam = A*(2*6.626e-27*(3e10)**2/(lam*1e-8)**5)/(np.exp(6.627e-27*3e10/(lam*1e-8*1.38e-16*T))-1)
    return Blam


# effective wavelengths of filters:
wle = {'S':2030,'D':2231,'P':2634,'u': 3560 , 'g': 4830 , 'r':6260 , 'i': 7670 , 'z': 9100 , 'U': 3600 , 'B': 4380 , 'V': 5450 , 'R': 6410 , 'I': 7980, 'J':12200, 'H':16300, 'K': 21900}


# Reference fluxes for converting magnitudes to flux
zp = {'S': 520 ,'D': 400 ,'P': 430 ,'u': 859.5 , 'g': 466.9 , 'r': 278.0 , 'i': 185.2 , 'z': 131.5 , 'U': 417.5 , 'B': 632 , 'V': 363.1 , 'R': 217.7 , 'I': 112.6, 'J':31.47, 'H':11.38, 'K':3.961}


bands = 'SDPugriz'

wl = []
fref = []

for i in bands:
    wl.append(wle[i])
    fref.append(zp[i]*1e-11) # because reference fluxes in units of 10^-11 erg/s/cm2/A

wl = np.array(wl)
fref = np.array(fref)

# note: can struggle to fit u/U sometimes, since a Supernova isn't a good BB in the UV!
mags = np.array([17.401,16.908,16.564,17.08,16.87,16.95,17.10,17.21])

z = 0.1

######## cosmocalc ############

# initialize constants
H0 = 72.                         # Hubble constant
WM = 0.27                        # Omega(matter)
WV = 1.0 - WM - 0.4165/(H0*H0)  # Omega(vacuum) or lambda

WR = 0.        # Omega(radiation)
WK = 0.        # Omega curvaturve = 1-Omega(total)
c = 299792.458 # velocity of light in km/sec
Tyr = 977.8    # coefficent for converting 1/H into Gyr
DTT = 0.0      # time from z to now in units of 1/H0
DCMR = 0.0     # comoving radial distance in units of c/H0
DA = 0.0       # angular size distance
DL = 0.0       # luminosity distance
DL_Mpc = 0.0
a = 1.0        # 1/(1+z), the scale factor of the Universe
az = 0.5       # 1/(1+z(object))

h = H0/100.
WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
WK = 1-WM-WR-WV
az = 1.0/(1+1.0*z)
n=1000         # number of points in integrals


for i in range(n):
    a = az+(1-az)*(i+0.5)/n
    adot = np.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
    DTT = DTT + 1./adot
    DCMR = DCMR + 1./(a*adot)

DTT = (1.-az)*DTT/n
DCMR = (1.-az)*DCMR/n

ratio = 1.00
x = np.sqrt(abs(WK))*DCMR
if x > 0.1:
    if WK > 0:
        ratio =  0.5*(np.exp(x)-np.exp(-x))/x
    else:
        ratio = np.sin(x)/x
else:
    y = x*x
    if WK < 0: y = -y
    ratio = 1. + y/6. + y*y/120.

DCMT = ratio*DCMR
DA = az*DCMT

DL = DA/(az*az)

DL_Mpc = (c/H0)*DL

#############################################

SN_distance = DL_Mpc*3.086e24 # convert Mpc to cm, since flux in erg/s/cm2/A

# convert mags to flux (and plot):

flux = 4*np.pi*SN_distance**2*fref*10**(-0.4*mags)
print 'Luminosities = ', flux
print "Flambda = ", fref*10**(-0.4*mags)
print "Fnu = ", (fref*10**(-0.4*mags))*(wl**2/(3e10))

plt.scatter(wl,flux)

# fit black body:

BBparams, covar = curve_fit(bbody,wl,flux,p0=(10000,4*np.pi*(1e15)**2))

# p0 is initial parameters to guess, to speed up finding fit.

T = BBparams[0]

Area = BBparams[1]

### remove the following comments if you want errors and chi2 ###
#resid = (flux - bbody(wl, *BBparams)) / len(flux)
#chisq = np.dot(resid, resid)
#print 'Chi squared = %.2f' % chisq

# BB temperature
print '\nBlackbody temperature = %.0f +\- %.0f K\n' % (T,np.sqrt(covar[0,0]))


# write the ouputs of the code in a .dat file
file = open("output.dat","w")

w = []
f = []


for wav in range(wl[0]-200,wl[-1]+200):
    file.write("%g\t%g\n" % (wav,bbody(wav,T,Area)))
    w.append(wav)
    f.append(bbody(wav,T,Area))

## option to give back the magnitudes estimated from the best BB fit. Here an example with g-band
# to be improved/automatised
if _magback == 'yes':
    fref1 = 466.9e-11
    flux1 = f[w.index(4830)]
    mag = log10(flux1/(4*np.pi*SN_distance**2*fref1))/(-0.4)
    print 'g= %0.3f' % (mag)

plt.plot(w,f,"r-")

w1 = []
f1 = []

for wav in range(2000,24000):
    file.write("%g\t%g\n" % (wav,bbody(wav,T,Area)))
    w1.append(wav)
    f1.append(bbody(wav,T,Area))

lum = max(integrate.cumtrapz(f1,w1))
print '#####################################################'
print 'Meaurements obtained with a BB from 2000 to 24000 Ang' # change input in line 206 if a different baseline is needed
print '#####################################################'
print '\nBolometric Luminosity = %.2E\n'  %lum

rr = sqrt(lum/(4*pi*(5.6704*10**-5)*T**4))
print '\nBlackbody radius = %.3E cm\n' %rr

peak = (2.8977729*10**7)/T
print '\nBlackbody peak = %.1f Angstrom\n' %peak
plt.show()


x = np.arange(1,10000)

y = bbody(x,T,Area)
plt.plot(x,y)
