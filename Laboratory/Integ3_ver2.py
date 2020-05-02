# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Laboratory

Date of creation: Thu Apr 30 01:01:10 2020

Version: 2 

Description:
        1. Defined the numerical integrator superclass
        2. Defined the Euler integrator subclass
        3. Definition of the point mass gravity function taken from ver1
        4. Defined the RK4 integrator subclass
        5. Output validation - through step size variation part of the assignment - validated
"""
###########################################
#### Importing the required modules  ######   
###########################################

import numpy as npMain
import spiceypy as spice
import sys

sys.path.append('../') # add everything that exists in MoMpy directory

from Modules.BasicAstrodynamics import convertKeplerToCartesian
from Modules.BasicAstrodynamics import convertCartesianToKepler

#################################################
## Definition: Numerical integrator superclass ##    
#################################################
class numericalIntegrator:
    """Definition of the superclass for a numerical integration operation"""
    
    # constructor
    def __init__(self,state0):        
        self.state0 = state0
        self._integrator = ""
        self.t0 = 0
        self.te = 0
        self.tStep = 0
    
    
    def setStateDerivativeFunction(self,sdf):
        self._sdf = sdf
                      
    def __basicIntegrationSetting(self,timeArray):
         
        # error handling
        assert(timeArray.size==3),"Incorrect size for the time array" 
        assert(te>=t0),"Final epoch of integration is smaller than intial epoch."
        
        # importing modules for the class definition
        import numpy as npNIC
        
        # extracting size of the state vector and time properties
        self.stateSize = self.state0.size
        self.t0 = timeArray[0]
        self.te = timeArray[1]
        self.tStep = timeArray[2]
    
        # Creating variables of interest
        self.time  = npNIC.arange(self.t0,self.te+self.tStep,self.tStep) # time vector
        Sol = npNIC.zeros(self.time.size*(1+self.stateSize)) 
        Sol.shape = [self.time.size,(1+self.stateSize)]
        
        # Initialization 
        Sol[0,0] = self.time[0] # saving intial epoch
        for index in range(1,self.stateSize):
            Sol[0,index] = self.state0[index-1] # saving initial state 
        
        return Sol
        
    # documenting the associated subclasses
    @staticmethod 
    def listAllSubclasses():
        subClist = ["Euler Integrator", "RK4 Integrator"]
        return subClist
        

###########################################
## Definition: Euler Integrator Subclass ##   
###########################################
class EulerIntegrator(numericalIntegrator):
    def integrate(self,timeArray):
        
        self._integrator = "Euler Integrator"
        
        Sol = super()._numericalIntegrator__basicIntegrationSetting(timeArray)
        
        # Using Euler integration formulation at each step
        for count in range (1, self.time.size):
            S = Sol[count-1,1:self.stateSize+1] # creates another array representing state at earlier epoch
        
            # Current epoch
            Sol[count,0] = self.time[count]
        
            # calling the state derivative function
            S_dot = self._sdf(S)
            
            S = S + S_dot * self.tStep #computing state at this epoch, using Euler formulation
            for index in range(1,self.stateSize):
                Sol[count,index] = S[index-1] # saving current state   
        
        return Sol
###########################################
### Definition: RK4 Integrator Class ###   
###########################################
class RK4Integrator(numericalIntegrator):
    def integrate(self,timeArray):
        self._integrator = "RK4 Integrator"
        
        Sol = super()._numericalIntegrator__basicIntegrationSetting(timeArray)
        
        # Using RK4 integration formulation at each step
        for count in range (1,self.time.size):
            S = Sol[count-1,1:self.stateSize+1] # creates another array representing state at earlier epoch
        
            # Current epoch
            Sol[count,0] = self.time[count]
        
            t1 = self.time[count-1] # last epoch
            K1 = self._sdf(S)
             
            t2 = t1+ step/2
            S2 = S + step* 0.5 * K1
            K2 = self._sdf(S2)
            
            t3 = t1 + step/2
            S3 = S + step* 0.5 * K2 
            K3 = self._sdf(S3)
     
            t4 = t1 + step
            S4 = S + step * K3 
            K4 = self._sdf(S4)
        
            phi = (K1 + 2*K2 + 2*K3 + K4)/6
            S = S + phi * step
            for index in range(1,self.stateSize):
                Sol[count,index] = S[index-1] # saving current state
        
        return Sol
            
    
    
###########################################
# Definition: Point Mass Gravity Function #   
###########################################
    
def pointMassGravity(S):
    "Function to compute point mass gravity acceleration and derivate function value for a state of the S/C"
    
    # error handling
    assert(S.size==6),"Incorrect size for the state argument of the pointMassGravity function"
    
    # importing required modules
    import math
    import numpy as np
    
    # extract the following from SPICE
    mu2 = 398600435436095.94 # mu for Earth
    
    # computing the acceleration
    r = math.sqrt(S[0]**2 + S[1]**2 + S[2]**2)
    S_dot = np.array([S[3], S[4], S[5], -mu2*S[0]/(r**3), -mu2*S[1]/(r**3), -mu2*S[2]/(r**3)])
        
    return S_dot

    
#####################################
########### User inputs #############    
#####################################
        
# Parameters definition
step = 1000  # [s]
t0 = 0
te= 30 * 86400 # [one month = 30 * sidereal day in seconds]
a = 42164173 # Semi-major axis for a GEO S/C [m]

KeplerC = npMain.array([a,0,0,0,0,0]) # initial state of GEO satellite, Kepler elements

tradeOff = 2 # 0 - No trade-off
             # 1 - Checking quality of results with the Euler integrator
             # 2-  Checking quality of results with the RK4 integrator

        
visualize = 0 # 0 - no visualizations
              # 1 - [UPDATE THIS]
              
######################################
######## Computations ################
######################################
              
# Importing the GM parameter from SPICE
spice.furnsh("../External_files/Spice_kernels/Load_Kernels.txt")              
muE = spice.bodvrd('EARTH','GM',1)
mu = muE[1][0]*1e9 # Gravitational parameter for Earth [m^3/s^2]
spice.kclear
             
# Converting from initial state in KE to equivalent cartesian coordinates
S0 = convertKeplerToCartesian(KeplerC,mu,3) # initial state in cartesian coordinates

if tradeOff == 0:
    # Integration using the Euler integrator
    Object1 = EulerIntegrator(S0) 
    Object1.setStateDerivativeFunction(pointMassGravity)
    #print(Object1.state0)
    result = Object1.integrate(npMain.array([t0,te,step]))
    #print(result)
    #print(Object1._integrator)
    out = numericalIntegrator.listAllSubclasses()
    
    # Integration using the RK4 integrator
    object2 = RK4Integrator(S0)
    object2.setStateDerivativeFunction(pointMassGravity)
    result2 = object2.integrate(npMain.array([t0,te,step]))

elif tradeOff == 1 or tradeOff == 2:
    # Setting values to be used
    referenceA = KeplerC[0]
    referenceE = KeplerC[1]
    steps = npMain.array([1000,100,10,1,0.1])
    
    #Storage variable
    qualityCheck = npMain.zeros(7*steps.size)
    qualityCheck.shape = [steps.size,7]
    qualityCheck[:,1] = referenceA
    qualityCheck[:,4] = referenceE    
    
    # Creating and setting up the required object
    if tradeOff==1:
        NIObject = EulerIntegrator(S0)
    else:
        NIObject = RK4Integrator(S0)
    NIObject.setStateDerivativeFunction(pointMassGravity)
    
    
    # Quality check through iterations
    for row in range(steps.size):
        # Extracting the step size and solving using the desired integrator
        step = steps[row]
        qualityCheck[row,0] = step
        Sol = NIObject.integrate(npMain.array([t0,te,step]))
            
        # Extracting the final state and converting it into the equivalent cartesian coordinates
        FinalS = Sol[-1,:]
        FinalS_carte = convertCartesianToKepler(FinalS[1:7],mu)
        
        # Quality-check
        qualityCheck[row,2] = FinalS_carte[0]
        qualityCheck[row,3] = (qualityCheck[row,2]-qualityCheck[row,1])*100/qualityCheck[row,1] # % error in semi-major axis
        qualityCheck[row,5] = FinalS_carte[1]
        qualityCheck[row,6] = (qualityCheck[row,5]-qualityCheck[row,4]) # absolute error in eccentricity
        print('Quality-check complete at step size of ',step, ' s.')



######################################
######## Visualization ###############
######################################

if visualize == 0:
    print('No visualizations are requested by the user.')
elif visualize == 1:
    print('nothing yet')
    # [EDIT THIS]
elif visualize ==2:
    print('nothing yet')
    # [EDIT THIS]  
else:
    print('Unrecognized input for the visualization option.')

