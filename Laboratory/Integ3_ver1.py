# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Laboratory

Date of creation: Sat Apr 25 01:10:29 2020

Version: v1 (Base: Integ1 excercise)

Description:
        1. Integ-3 assignment to solve problem of propagation S/C state under point mass gravity
        2. Analytical calculation function from the base is dropped
        3. Function for the Euler integrator edited to call the derivative function
        4. State derivative function, for point mass gravity of Earth is defined
        5. Function for the RK4 integrator edited to call the derivative function
        6. Orbit visualizations for results obtained with both integrators
        7. Output validation - TBD
"""
###########################################
#### Importing the required modules  ######   
###########################################

import numpy as npMain
import matplotlib.pyplot as plt
import sys
import spiceypy as spice

sys.path.append('../') # one way to add everything that exists in MoMpy directory

from Modules.BasicAstrodynamics import convertCartesianToKepler
from Modules.BasicAstrodynamics import convertKeplerToCartesian

###########################################
## Definition: Euler Integrator Function ##   
###########################################

def Euler(step,t0,te,S):
    # importing the required modules
    import numpy as np
      
    # extracting size of the state vector
    stateSize = S.size
    
    # Creating variables of interest
    time  = np.arange(t0,te+step,step) # time vector
    EulerSol = np.zeros(time.size*(1+stateSize)) 
    EulerSol.shape = [time.size,(1+stateSize)]
    
    # Initialization 
    EulerSol[0,0] = time[0] # saving intial epoch
    for index in range(1,stateSize):
        EulerSol[0,index] = S[index-1] # saving initial state 
         
    # Using Euler integration formulation at each step
    for count in range (1, time.size):
        S = EulerSol[count-1,1:stateSize+1] # creates another array representing state at earlier epoch
        
        # Current epoch
        EulerSol[count,0] = time[count]
        
        # calling the state derivative function
        S_dot = pointMassGravity(S)
            
        S = S + S_dot * step #computing state at this epoch, using Euler formulation
        for index in range(1,stateSize):
            EulerSol[count,index] = S[index-1] # saving current state         
        
    return EulerSol

###########################################
### Definition: RK4 Integrator Function ###   
###########################################
    
def RK4(step,t0,te,S):
    # importing the required modules
    import numpy as np
    
    # extracting the size of state vector
    stateSize = S.size
    
    # Creating variables of interest
    time = np.arange(t0,te+step,step) # time vector
    RK4Sol = np.zeros(time.size*(1+stateSize)) 
    RK4Sol.shape = [time.size,(1+stateSize)] # storage for solution

    # Initialization 
    RK4Sol[0,0] = time[0] # saving intial epoch
    for index in range(1,stateSize):
        RK4Sol[0,index] = S[index-1] # saving initial state 
    
    # Using RK4 integration formulation at each step
    for count in range (1,time.size):
        S = RK4Sol[count-1,1:stateSize+1] # creates another array representing state at earlier epoch
        
        # Current epoch
        RK4Sol[count,0] = time[count]
        
        t1 = time[count-1] # last epoch
        K1 = pointMassGravity(S)
             
        t2 = t1+ step/2
        S2 = S + step* 0.5 * K1
        K2 = pointMassGravity(S2)
            
        t3 = t1 + step/2
        S3 = S + step* 0.5 * K2 
        K3 = pointMassGravity(S3)
     
        t4 = t1 + step
        S4 = S + step * K3 
        K4 = pointMassGravity(S4)
        
        phi = (K1 + 2*K2 + 2*K3 + K4)/6
        S = S + phi * step
        for index in range(1,stateSize):
            RK4Sol[count,index] = S[index-1] # saving current state
    
    return RK4Sol

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
step = 10  # [s]
t0 = 0
te= 30 * 86400 # [one month = 30 * sidereal day in seconds]

KeplerC = npMain.array([42164173,0,0,0,0,0]) # initial state of GEO satellite, Kepler elements

visualize = 2 # 0 - no visualizations
              # 1 - plot the S/C orbit over a month, XY plane, Euler integrator
              # 2 - plot the S/C orbit over a month, XY plane, RK4 integrator
           
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

# Solving the problem using Euler integrator
EulerSol = Euler(step,t0,te,S0) 

# Solving the problem using RK4 integrator
RK4Sol = RK4(step,t0,te,S0) # results validated at multiple step-sizes

######################################
######## Visualization ###############
######################################

if visualize == 0:
    print('No visualizations are requested by the user.')
elif visualize == 1:
    # Plotting the orbit of S/C in X-Y plane, when Euler integrator is used
    plt.plot(0.001 * EulerSol[:,1],0.001 * EulerSol[:,2])
    plt.title('Orbit of the S/C over a period of one month, when the EOM is integrated using the Euler integrator')
    plt.xlabel('X [km]')
    plt.ylabel('Y [km]')
    plt.grid(True)
elif visualize ==2:
    # Plotting the orbit of S/C in X-Y plane, when RK4 integrator is used
    plt.plot(0.001 * RK4Sol[:,1],0.001 * RK4Sol[:,2])
    plt.title('Orbit of the S/C over a period of one month, when the EOM is integrated using the RK4 integrator')
    plt.xlabel('X [km]')
    plt.ylabel('Y [km]')
    plt.grid(True)    
else:
    print('Unrecognized input for the visualization option.')

