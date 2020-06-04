# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Modules

Date of creation: Wed Jun  3 16:09:06 2020

Version: Base 

Description:
        1. Module to contain useful astrodynamic problems
        2. LEO to GEO two-burn transfer problem
        3. 
        4. 
        5. 
        6. 
        7. 
"""
###############################################################
#### Importing the required modules - for local testing  ######   
###############################################################

import numpy as np

###########################################
## Definition: Two-burn trasfer problem ###   
###########################################
class TwoBurnTransfer:
    """ Defining two burn transfer problem from/to LEO to/from GEO """
    
    # constructor - assigning and checking the problem type
    def __init__(self,problemtype):
        self.problemType = problemtype 
        self.__paraSet = False
        
        # error handling
        assert(self.problemType==1),"Currently, only problem type 1 is supported. Check the descriptions of supported problems using 'stateProblemTypes' method of the class."

        # setting up the number of required parameters, to be provided by user
        if self.problemType==1:
            self.__requiredPara = 2
    
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
            
            
        # Defining/setting parameters that are specific to the problem-type
        if self.problemType==1:
            # Importing required modules
            import math
            
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
        
            self.__paraSet = True # Parameters are not completely set and other methods can be used
            
        else:
            print('Problem type is not supported.')
        
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
            
        else:
            print('Problem type not supported.')
            
        return DeltaVs
    
    # Method to compute the required derivative value
    def computeDerivative(self,variablesSet):    
        # Error handling
        assert(self.__paraSet == True),"First set the parameters of the object using the method - 'setProblemParameters'."
        
        # importing the required modules
        import math
        
        if self.problemType==1:
            # computations
            i1 = variablesSet * math.pi/180
            i2 = (self.__iLEO - variablesSet) * math.pi/180
            
            A = self.__VLEO * self.__Vp * math.sin(i1) / math.sqrt(self.__VLEO**2 + self.__Vp**2 - 2*self.__VLEO*self.__Vp*math.cos(i1))
            B = self.__VGEO * self.__Va * math.sin(i2) / math.sqrt(self.__VGEO**2 + self.__Va**2 - 2*self.__VGEO*self.__Va*math.cos(i2))
            
            return (A-B)
    
    # Documenting problem types
    @staticmethod 
    def stateProblemTypes():
        problems = ["Problem type 1: Circular LEO to GEO"]
        return problems

              
######################################
######## Computations ################
######################################


######################################
######## Local testing ###############
######################################
object1 = TwoBurnTransfer(1)
#print(object1.problemType)

para = np.array([185,28.5])
object1.setProblemParameters(para)
result = object1.computeDeltaV(10)
result2 = object1.computeDerivative(2.165)


#temp2 = TwoBurnTransfer.stateProblemTypes()



