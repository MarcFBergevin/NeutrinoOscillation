#
# Program to obtain the neutrino spectrum as a function of stand-off L
#

#
# Imports
#
import math as math
import numpy as n
import ROOT as ROOT
from ROOT import TCanvas, TPad, TFile, TPaveText
from ROOT import gBenchmark, gStyle, gROOT,TColor, TChain

#
# Classes definition
#
class NeutrinoOscillation:
    """ Class to display the neutrino oscillation at a specific stand-off """

    #
    # Class Members
    #
    rate        = None
    s12			= 0.0
    s23			= 0.0
    s13			= 0.0
    dm12 		= 0.0
    dm23		= 0.0
    dm13 		= 0.0
    t12 		= 0.0
    t23 		= 0.0
    t13 		= 0.0
    delta_31 	= 0.0
    delta_21 	= 0.0
    delta_32 	= 0.0

    #
    # Class methods
    #
    def __init__(self, medium,mass,power,standoff,
                    s12_flag=None,s23_flag=None,s13_flag=None,
                    dm12_flag=None,dm23_flag=None,hierarchy_flag=0):
        """ Initialise all the default paramters"""
        """ Values taken from http://pdg.lbl.gov/2013/reviews/rpp2013-rev-neutrino-mixing.pdf """
    
        if dm12_flag == None:
            dm12 		= 7.54e-5               # (+0.26/-0.22)
        if s12_flag == None:
            s12			= math.sqrt(0.307)      # (+0.018/-0.018)
        
        if   hierarchy_flag == 0:
            dm13 		= 2.43e-3 + dm12/2.0        # (+0.07/-0.11)
            dm23        = 2.43e-3 - dm12/2.0
            if s23_flag == None:
                s23			= math.sqrt(0.386)  # (+0.024/-0.024)
            if s13_flag == None:
                s13			= math.sqrt(0.0241) # (+0.0025/-0.0025)
            print '\nUsing normal hierarchy'

        elif hierarchy_flag == 1:
            dm13 		= 2.42e-3 + dm12/2.0        # (+0.07/-0.11)
            dm23        = 2.42e-3 - dm12/2.0
            if s23_flag == None:
                s23			= math.sqrt(0.392)  # (+0.024/-0.024)
            if s13_flag == None:
                s13			= math.sqrt(0.0244) # (+0.0023/-0.0025)
            print '\nUsing inverted hierarchy'

        delta_31 	= 1.27*dm13
        delta_21 	= 1.27*dm12
        delta_32 	= 1.27*dm23
        t12 		= math.asin(s12)
        t23 		= math.asin(s23)
        t13 		= math.asin(s13)

        print ' %s_12 = %0.4e \t\t %s_12 = %4.2f%s \t\t sin^2(2%s_12) = %0.3f \t\t tan^2(%s_12)  = %0.3f' % (u"\u03B8",
                                                                                                             t12,u"\u03B8",
                                                                                                             t12*180./math.pi,u"\xb0",
                                                                                                             u"\u03B8",
                                                                                                             math.sin(2.0*t12)*math.sin(2.0*t12),
                                                                                                             u"\u03B8",
                                                                                                             math.tan(t12)*math.tan(t12))

        print ' %s_23 = %0.4e \t\t %s_23 = %4.2f%s \t\t sin^2(2%s_23) = %0.3f' % (u"\u03B8",t23,
                                                                                  u"\u03B8",
                                                                                  t23*180./math.pi,
                                                                                  u"\xb0",
                                                                                  u"\u03B8",
                                                                                  math.sin(2.0*t23)*math.sin(2.0*t23))

        print ' %s_13 = %0.4e \t\t %s_13 =  %4.2f%s \t\t sin^2(2%s_13) = %0.3f' % (u"\u03B8",
                                                                                   t13,
                                                                                   u"\u03B8",
                                                                                   t13*180./math.pi,
                                                                                   u"\xb0",
                                                                                   u"\u03B8",
                                                                                   math.sin(2.0*t13)*math.sin(2.0*t13))

        print ' %s_12 = %0.4e eV**2' % (u"\u0394",dm12)
        print ' %s_23 = %0.4e eV**2' % (u"\u0394",dm23)
        print ' %s_13 = %0.4e eV**2\n' % (u"\u0394",dm13)

        GWth       = 2e20              # neutrino per GWth
        av_num     = 6.02214129e23     # mol^-1 2 H
        length     = standoff * 1e5    # cm/km
        xsect      = 9.54e-44          #
        time       = 60.0*60.0*24.0    # s/day
        if medium == 1:
            cubicCM    = 1e+6          # m^3 in cm^3
            dens       = 1.0           # g/cm3
            mol_weight = 2./18.02      # 2 mols of hydrogen in water
        elif medium ==0:
            cubicCM    = 1e+6          # m^3 in cm^3
            dens       = 0.804         # g/cm3 from paper
            mol_weight = 48./352.      #  mols of hydrogen in medium (52./352.not correct for DC)
        else:
            raise ValueError("You have chosen a non-existing material, Goodbye".format(medium))
        
        #molecular_weight = 52./352.# Pseudcomine

        protons    = av_num*mass*cubicCM*dens*mol_weight
        
        rate   = power * GWth * xsect * protons  * time / length**2/ (4.0*math.pi)

        print 'Detector and reactor characteristics:'
        print ' Detector mass             : %4.3f' %(mass)
        print '          proton-target    : %4.3e' %(protons)
        print ' Reactor  power (GWth)     : %4.3f' %(power)
        print '          distance (km)    : %4.3f' %(standoff)
        print ' Neutrino rate in detector : %4.3f per day (pre-efficiency)' %(rate)


    def FindRate(self,medium,mass,power,standoff):
    

    
        return self.rate




#
# If run as script
#
if __name__ == "__main__":
    
    personalDetector    = raw_input('Would you like to define your own detector (yes/no)?:')
    
    if personalDetector == 'yes':
        # Define your very own detector
        detectorMedium      = int(raw_input('What detector medium (liquid scintillator=0,water=1)?:\n'))
        detectorMass        = int(raw_input('What detector mass (ton)?:\n'))
        reactorPower        = float(raw_input('What is the reactor power (GWth)?:\n'))
        reactorStandoff     = float(raw_input('What is the reactor to detector distance (km)?:\n'))
    else:
        # Define the Double Chooz detector as default
        detectorMedium      = 0
        detectorMass        = 10.3
        reactorPower        = 8.5
        reactorStandoff     = 1.05

    nuOsc = NeutrinoOscillation(detectorMedium,detectorMass,reactorPower,reactorStandoff,hierarchy_flag=0)
    # One can change any paramters manually in the following way
    #nuOsc = NeutrinoOscillation(detectorMedium,detectorMass,reactorPower,reactorStandoff,s13=3)
