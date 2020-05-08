# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Laboratory

Date of creation: Tue May  5 02:46:06 2020

Version: 3 

Description:
        1. Definition of NI classes imported from the Integrators module
        2. EOM are imported using the Propagators modules and passed as the state derivative function
        3. Output validation -TBD
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

from Modules.Integrators import * # not a recommendeed method?
from Modules.Propagators import Propagator 
    
#####################################
########### User inputs #############    
#####################################
        
# Parameters definition
step = 1000  # [s]
t0 = 0
te= 30 * 86400 # [one month = 30 * sidereal day in seconds]
a = 42164173 # Semi-major axis for a GEO S/C [m]

# Defining the propagtor settings 
prop = Propagator("Cowell")
prop.setCentralBody("Earth")
EOM = prop.pointMassGravity

KeplerC = npMain.array([a,0,0,0,0,0]) # initial state of GEO satellite, Kepler elements

tradeOff = 1 # 0 - No trade-off
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
    Object1.setStateDerivativeFunction(EOM)
    result = Object1.integrate(npMain.array([t0,te,step]))
    out = numericalIntegrator.listAllSubclasses()
    
    # Integration using the RK4 integrator
    object2 = RK4Integrator(S0)
    object2.setStateDerivativeFunction(EOM)
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
    NIObject.setStateDerivativeFunction(EOM)
    
    
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

