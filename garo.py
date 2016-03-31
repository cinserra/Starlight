#! /usr/bin/python

# Created by Cosimo Inserra
# GitHub: https://github.com/cinserra
# @COSMO_83

################### for the help ##################
from optparse import OptionParser
from numpy import *

description = "Script to evaluate the effective gain and readout for multiple, combined images"
usage = "%prog "
if __name__ == "__main__":
    parser = OptionParser(usage=usage, description=description,)
    parser.add_option("-v", "--verbose",dest="verbose",\
                  action="store_true",default=False,
                  help='Print tasks description')
    parser.add_option("-t", "--telescope",dest="telescope", action="store",default='', type="str",
                  help='Select the telescope (NTT-E/NTT-S/LT)')
    parser.add_option("-f", "--frames",dest="nframe", action="store",type="float" ,default=None,
                  help='Number of frames')
    parser.add_option("-c", "--combine",dest="combine", action="store",default='', type="str",
                  help='Select the combination method (media/average/sum)')
    option,args = parser.parse_args()



_telescope = option.telescope

if _telescope == 'NTT-E': #EFOSC camera
	gain = 1.18
	ron = 11.6
elif _telescope == 'NTT-S': #SOFI camera
	gain = 5.4
	ron = 12.0
elif _telescope == 'LT': #IO camera
	gain = 1.62
	ron = 8.0

number = option.nframe
_combine = option.combine
#info from DAOPHOT primer
if _combine == 'median':
	effgain = 2. * number * gain / 3.
	effron = sqrt(2./3. * number) * ron
elif _combine == 'average':
	effgain = float(number) * float(gain)
	effron = sqrt(number) * ron
elif _combine == 'sum':
	effgain = float(gain)
	effron = sqrt(number) * ron

print '####################################'
print ''
print _telescope, 'with %0.0f images' % (number)
print 'Effective gain = %0.2f' % (effgain)
print 'Effective readout noise = %0.2f' % (effron)
print ''
print '####################################'
