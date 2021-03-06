# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Modules

Date of creation: Tue Jun 30 10:38:29 2020

Version: 1 

Description:
        1. Module to contain useful astrodynamic problems
        2. W.r.t. base class, removed the local testing of transfer orbot problems of type 1 and 2
"""


###########################################
## Definition: Two-burn trasfer problem ###   
###########################################
class TransferOrbits:
    """ Defining two/three burn transfer problem from/to LEO to/from GEO/Higher altitude orbit """
    
    ###########################
    ###########################  
    
    # constructor - assigning and checking the problem type
    def __init__(self,problemtype):
        self.problemType = problemtype 
        self.__paraSet = False
        
        # error handling
        assert(self.problemType==1 or self.problemType==2),"Currently, only problem type 1 and 2 are supported. Check the descriptions of supported problems using 'stateProblemTypes' method of the class."

        # setting up the number of required parameters, to be provided by user
        if self.problemType==1:
            self.__requiredPara = 2
        elif self.problemType == 2:
            self.__requiredPara = 4
    
    ###########################
    ###########################  
    
    # Setting up the parameters of the problem
    def setProblemParameters(self,ParaSet):   
        self.__Paraset = ParaSet
        
        # Error-handling
        assert(self.__Paraset.size==self.__requiredPara),"Incorrect array-size for the parameter matrix."
        
        
        # Extracting required parameters using SPICE kernels if possible
        try:
            import spiceypy as spice
            import sys

            sys.path.append('../') # add everything that exists in the MoMPy directory
            
            spice.furnsh('../External_files/Spice_kernels/Load_Kernels3.txt')
            muE = spice.bodvrd('EARTH','GM',1)
            mu = muE[1][0] # mu for Earth [km^3/s^2]
            RE = spice.bodvrd('EARTH','RADII',3)
            Re= RE[1][0] # [km], Average radius of Earth
            spice.kclear    
        except:
            mu = 398600.435436095925979316234588623046875 # manually entering in case the import fails
            Re = 6378.136599999999816645868122577667236328125 # manually entering in case the import fails
        finally:
            hGEO = 35768 # Altitude of GEO [km]
            
        # Importing required modules
        import math
        
        # Defining/setting parameters that are specific to the problem-type
        if self.problemType==1:
            # Extraction
            hLEO = self.__Paraset[0] # Altitude of circular LEO [km]
            self.__iLEO = self.__Paraset[1] # Inclination of LEO [deg]
            assert( hLEO < hGEO ), "LEO altitude can not be greater than GEO altitude. Check inputs."
            
            # Computation
            ra = Re + hGEO # apoapsis of the GTO [km]
            rp = Re + hLEO # periapsis of the GTO [km]
            self.__VLEO = math.sqrt(mu/rp)  # [km/s]
            self.__VGEO = math.sqrt(mu/ra) # [km/s]          
            e = (ra-rp)/(ra+rp) # eccentricity of the GTO
            a = (ra +rp)/2 # semi-major axis of the GTO
            self.__Va = math.sqrt(mu * (1-e)/(a * (1+e))) # Apoapsis velocity in the GTO [km/s] 
            self.__Vp = math.sqrt(mu * (1+e)/(a * (1-e))) # Periapsis velocity in the GTO [km/s]
        
            self.__paraSet = True # Parameters are now completely set and other methods can be used
        elif self.problemType == 2:
            # Extraction
            self.__r1LEO = Re + self.__Paraset[0] # Semi-major axis of the circular LEO [km]
            self.__iLEO = self.__Paraset[1] # Inclination of the circular LEO [deg]
            self.__r2HEO = Re + self.__Paraset[2] # Semi-major axis of the high altitude final orbit [km]
            self.__raInterMedLimit = Re + self.__Paraset[3] # Limiting value of the apoapsis distance of the intermediate transfer orbit [km]
            assert( self.__r2HEO < self.__raInterMedLimit ), "Limiting value of the apoapsis distance of the intermediate transfer orbit has to be greater than semi-major axis of the final circular-equatorial orbit. Check inputs."
            
            # Computation
            self.__VC1 = math.sqrt(mu/self.__r1LEO) # Velocity of the circular LEO [km/s]
            self.__VC2 = math.sqrt(mu/self.__r2HEO) # Velocity of the circular final orbit [km/s]            
            
            self.__paraSet = True # Parameters are now completely set and other methods can be used
        else:
            print('Problem type is not supported.')
            
    ###########################
    ###########################        
        
    # Method to compute Delta Vs for the given set of user inputs
    def computeDeltaV(self,variablesSet):
        
        # Error handling
        assert(self.__paraSet == True),"First set the parameters of the object using the method - 'setProblemParameters'."
        
        # importing the required modules
        import math
        import numpy as np
        
        if self.problemType==1:
            # storage for the results
            DeltaVs = np.zeros(2)
            
            # computations
            i1 = variablesSet * math.pi/180
            i2 = (self.__iLEO - variablesSet) * math.pi/180
            
            DeltaVs[0] = math.sqrt(self.__VLEO**2 + self.__Vp**2 - 2*self.__VLEO*self.__Vp*math.cos(i1)) # DeltaV for first burn
            DeltaVs[1] = math.sqrt(self.__VGEO**2 + self.__Va**2 - 2*self.__VGEO*self.__Va*math.cos(i2)) # DeltaV for second burn
        
        elif self.problemType==2:
            # storage for the results
            DeltaVs = np.zeros(3)
            
            # Error-handling
            assert(variablesSet[1]>0),"Check the second input variable - only positive values are accepted."
            
            # computations
            i2 = variablesSet[0] * math.pi/180 # Inclination change at the second burn
            i3 = (self.__iLEO - variablesSet[0]) * math.pi/180 # Inclination change at the third burn
            ra = variablesSet[1] * (self.__raInterMedLimit-self.__r2HEO) # Apoapsis altitude of the intermediate transfer orbit [km]
                        
            et1 = (ra-self.__r1LEO)/(ra + self.__r1LEO) # Eccentricity of the first transfer leg 
            et2 = (ra-self.__r2HEO)/(ra + self.__r2HEO) # Eccentricity of the second transfer leg
            VCa = self.__VC1 * math.sqrt(self.__r1LEO/ra) # Velocity of the circular orbit with radius same as the apoapsis distance of the intermediate transfer orbit
            
            # First Delta-V required, to go from circular LEO to first leg of the intermediate elliptical transfer orbit, orbital plane is not changed 
            DeltaVs[0] = self.__VC1 * (math.sqrt(1+et1)-1)
            # Second Delta-V required, at the apoapsis of the first leg of the intermediate ellitptical transfer orbit, to put S/C on the second leg, an inclination change is achieved
            DeltaVs[1] = VCa * math.sqrt( (1-et1) + (1-et2) -2 * math.sqrt((1-et1)*(1-et2)) * math.cos(i2) )
            # Third Delta-V required, at the periapsis of the second leg of the intermediate ellitptical transfer orbit, to put S/C in the final circular orbit, an inclination change is achieved
#            DeltaVs[2] = self.__VC2 * (math.sqrt(1+et2)-1) #Original formula in wakker when there is no inclination change at this maneuver
            DeltaVs[2] = self.__VC2 * math.sqrt( (1+et2) + 1 - 2 *math.sqrt(1+et2) * math.cos(i3))
        else:
            print('Problem type not supported.')
            
        return DeltaVs
    
    ###########################
    ###########################  
    
    # Method to compute the required derivative value
    def computeDerivative(self,variablesSet):    
        # Error handling
        assert(self.__paraSet == True),"First set the parameters of the object using the method - 'setProblemParameters'."
        
        # importing the required modules
        import math
        import numpy
        
        if self.problemType==1:
            # computations
            i1 = variablesSet * math.pi/180
            i2 = (self.__iLEO - variablesSet) * math.pi/180
            
            A = self.__VLEO * self.__Vp * math.sin(i1) / math.sqrt(self.__VLEO**2 + self.__Vp**2 - 2*self.__VLEO*self.__Vp*math.cos(i1))
            B = self.__VGEO * self.__Va * math.sin(i2) / math.sqrt(self.__VGEO**2 + self.__Va**2 - 2*self.__VGEO*self.__Va*math.cos(i2))
            
            return (A-B)
        elif self.problemType==2:
            # Error-handling
            assert(variablesSet[1]>0),"Check the second input variable - only positive values are accepted."
            
            # creating a storage variable
            derivatives = numpy.zeros(6)
            derivatives.shape = [3,2]
            
            #Extraction
            i2 = variablesSet[0] * math.pi/180 # Inclination change at the second burn
            i3 = (self.__iLEO - variablesSet[0]) * math.pi/180 # Inclination change at the third burn
            f = variablesSet[1] # Fraction determining the apoapsis altitute of the intermediate transfer orbit
            ra = f * (self.__raInterMedLimit-self.__r2HEO) # Apoapsis altitude of the intermediate transfer orbit [km]
            
            # Computations
            et1 = (ra-self.__r1LEO)/(ra + self.__r1LEO) # Eccentricity of the first transfer leg 
            et2 = (ra-self.__r2HEO)/(ra + self.__r2HEO) # Eccentricity of the second transfer leg
            
            # Defining essential functions
            f1 = 1 + et1
            f2 = 2 - f1
            g2 = 1 - et2
            g1 = 2 - g2
            h1 = self.__VC1 * math.sqrt(self.__r1LEO/ra)
            k1 = math.sqrt( f2 + g2 - 2 * math.sqrt(f2*g2) * math.cos(i2))
            k2 = math.sqrt( g1 + 1 - 2 * math.sqrt(g1) * math.cos(i3))
            
            # Defining the derivatives of above essential functions
            df1df = (2*self.__r1LEO*(self.__raInterMedLimit-self.__r2HEO))/(f*(self.__raInterMedLimit-self.__r2HEO)+self.__r1LEO)**2
            dh1df = -0.5 * self.__VC1 * math.sqrt(self.__r1LEO/(self.__raInterMedLimit-self.__r2HEO)) * math.pow(f,-1.5)
            df2df = -df1df
            dg2df = (-2 *self.__r2HEO*(self.__raInterMedLimit-self.__r2HEO)) / (f*(self.__raInterMedLimit-self.__r2HEO)+self.__r2HEO)**2 
            dg1df = - dg2df 
            dk1df = 0.5 * math.pow(k1,-3) * (df2df * dg2df - math.cos(i2) * math.pow(f2*g2,-1.5) * (f2*dg2df+df2df*g2) )
            
            # Derivatives with respect to first inclination change- i2
            # Derivative of DeltaV1
            derivatives[0,0]= 0
            # Derivative of DeltaV2
            derivatives[1,0]= h1 * math.pow(k1,-3) * math.sqrt(f2*g2) * math.sin(i2)
            # Derivative of DeltaV3
            derivatives[2,0]= - self.__VC2 * math.pow(k2,-3) * math.sqrt(g1) * math.sin(i3)
            
            # Derivatives with respect to fraction f
            # Derivative of DeltaV1
            derivatives[0,1]= 0.5 * self.__VC1 * df1df/math.sqrt(f1) 
            # Derivative of DeltaV2
            derivatives[1,1]= dh1df * k1 + h1 * dk1df  
            # Derivative of DeltaV3
            derivatives[2,1]= 0.5 * self.__VC2 * math.pow(k2,-3) * ( g1 - math.pow(g1,-1.5) * dg1df * math.cos(i3) ) 
            
            return derivatives
        
        
    # Documenting problem types
    @staticmethod 
    def stateProblemTypes():
        problems = ["Problem type 1: Circular LEO to GEO, two-burn transfer", "Problem type 2: Circular LEO to Higher altitude circular-equatorial orbit, three-burn transfer"]
        return problems

              
