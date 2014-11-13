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
    # Class Members.  Add
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
    # Define a ROOT object
    #
    U235S  = ROOT.TF1("U235S","[0]*exp([1]-[2]*x-[3]*x*x)",1.806,20.00)
    Pu239S = ROOT.TF1("Pu239S","[0]*exp([1]-[2]*x-[3]*x*x)",1.806,14.000)
    U238S  = ROOT.TF1("U238S","[0]*exp([1]-[2]*x-[3]*x*x)",1.806,14.000)
    Pu241S = ROOT.TF1("Pu241S","[0]*exp([1]-[2]*x-[3]*x*x)",1.806,14.000)
    IBDEnergyTotal = ROOT.TF1
    
    nu_U235S   = ROOT.TF1
    nu_Pu239S  = ROOT.TF1
    nu_U238S   = ROOT.TF1
    nu_Pu241S  = ROOT.TF1

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
            mol_weight = 48./355.      #  mols of hydrogen in medium (52./352.not correct for DC)
        else:
            raise ValueError("You have chosen a non-existing material, Goodbye".format(medium))
        
        #molecular_weight = 52./352.# Pseudcomine

        protons    = av_num*mass*cubicCM*dens*mol_weight
        
        self.rate   = power * GWth * xsect * protons  * time / length**2/ (4.0*math.pi)

        print 'Detector and reactor characteristics:'
        print ' Detector mass             : %4.3f' %(mass)
        print '          proton-target    : %4.3e' %(protons)
        print ' Reactor  power (GWth)     : %4.3f' %(power)
        print '          distance (km)    : %4.3f' %(standoff)
        print ' Neutrino rate in detector : %4.3f per day (pre-efficiency)' %(self.rate)

        mNeutron        = 939.565378
        mProton         = 938.27
        mElectron       = 0.51099891
        delta           = mNeutron - mProton
        A               = 0.5
        B               = mNeutron*mNeutron
        C               = 4.0*mProton
        D               = delta+(delta*delta-mElectron*mElectron)/(2*mProton)
        E               = mNeutron

        alpha_235U      = 0.496
        alpha_239P      = 0.351
        alpha_238U      = 0.087
        alpha_241P      = 0.066

        beta_235U      = 0.870
        beta_239P      = 0.896
        beta_238U      = 0.976
        beta_241P      = 0.793
        
        gamma_235U      = 0.160
        gamma_239P      = 0.239
        gamma_238U      = 0.162
        gamma_241P      = 0.080

        kappa_235U      = 0.091
        kappa_239P      = 0.0981
        kappa_238U      = 0.079
        kappa_241P      = 0.1085


        Ef_235U         = 201.92
        Ef_239P         = 209.99
        Ef_238U         = 205.52
        Ef_241P         = 213.60
    
        Ef = alpha_235U*Ef_235U + alpha_239P*Ef_239P + alpha_238U*Ef_238U + alpha_241P*Ef_241P
        
        self.U235S.SetParameter(0,alpha_235U)
        self.U235S.SetParameter(1,beta_235U)
        self.U235S.SetParameter(2,gamma_235U)
        self.U235S.SetParameter(3,kappa_235U)

        self.Pu239S.SetParameter(0,alpha_239P)
        self.Pu239S.SetParameter(1,beta_239P)
        self.Pu239S.SetParameter(2,gamma_239P)
        self.Pu239S.SetParameter(3,kappa_239P)

        self.U238S.SetParameter(0,alpha_238U)
        self.U238S.SetParameter(1,beta_238U)
        self.U238S.SetParameter(2,gamma_238U)
        self.U238S.SetParameter(3,kappa_238U)

        self.Pu241S.SetParameter(0,alpha_241P)
        self.Pu241S.SetParameter(1,beta_241P)
        self.Pu241S.SetParameter(2,gamma_241P)
        self.Pu241S.SetParameter(3,kappa_241P)

        IBDXsect = ROOT.TF1("IBDXsect","[0]*(sqrt([1]-[2]*([3]-x))-[4])",1.806,14.000)
        IBDXsect.SetParameter(0,A)
        IBDXsect.SetParameter(1,B)
        IBDXsect.SetParameter(2,C)
        IBDXsect.SetParameter(3,D)
        IBDXsect.SetParameter(4,E)
        Xsect = ROOT.TF1("Xsect","IBDXsect*sqrt(IBDXsect**2-[0]**2)",1.806,14.000)
                                        
        self.IBDEnergy = ROOT.TF1("IBDEnergy","(U235S + Pu239S+U238S+Pu241S)*Xsect",1.806,14.000)
        self.nu_U235S   = ROOT.TF1("nu_U235S",  "(U235S)*Xsect",1.806,14.000)
        self.nu_Pu239S  = ROOT.TF1("nu_Pu239S","(Pu239S)*Xsect",1.806,14.000)
        self.nu_U238S   = ROOT.TF1("nu_U238S",  "(U238S)*Xsect",1.806,14.000)
        self.nu_Pu241S  = ROOT.TF1("nu_Pu241S","(Pu241S)*Xsect",1.806,14.000)
    
        self.IBDEnergy.SetNpx(10000)
        self.nu_U235S.SetNpx(10000)
        self.nu_Pu239S.SetNpx(10000)
        self.nu_U238S.SetNpx(10000)
        self.nu_Pu241S.SetNpx(10000)

        self.IBDEnergy.SetTitle("All isotopes")
        self.nu_U235S.SetTitle("U235 isotope")
        self.nu_Pu239S.SetTitle("Pu239 isotope")
        self.nu_U238S.SetTitle("U238 isotope")
        self.nu_Pu241S.SetTitle("Pu241 isotope")
    
            
    def ApplyBurnup(self,days):
        """ Add the burn-up behavior applied """
        
        alpha_235U_bu = alpha_235U 
        alpha_239P_bu = alpha_239P
        alpha_238U_bu = alpha_238U
        alpha_241P_bu = alpha_241P
        
        self.U235S.SetParameter(0,  alpha_235U_bu)
        self.Pu239S.SetParameter(0, alpha_239P_bu)
        self.U238S.SetParameter(0,  alpha_238U_bu)
        self.Pu241S.SetParameter(0, alpha_241P_bu)
    
       

    def FindRate(self):
    
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



    nuOsc.nu_U235S.SetLineColor(ROOT.kMagenta)
    nuOsc.nu_Pu239S.SetLineColor(ROOT.kOrange+8)
    nuOsc.nu_U238S.SetLineColor(ROOT.kBlue-6)
    nuOsc.nu_Pu241S.SetLineColor(ROOT.kGreen+3)
    print '%s ' %(nuOsc.FindRate())
    nuOsc.IBDEnergy.Draw()
    nuOsc.nu_U235S.Draw("same")
    nuOsc.nu_Pu239S.Draw("same")
    nuOsc.nu_U238S.Draw("same")
    nuOsc.nu_Pu241S.Draw("same")