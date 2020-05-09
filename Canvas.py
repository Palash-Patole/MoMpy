# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Laboratory
Date of creation: Fri Nov  8 20:57:12 2019

version: Base

Description:
    1. Generic canvas to test different modules
    2. Tested - coordinate coversations
    3. 3 options of computing the Earth-repeat orbit solutions
    4. 
"""

from Modules.BasicAstrodynamics import convertCartesianToKepler
from Modules.BasicAstrodynamics import convertKeplerToCartesian

from Modules.OrbitDesign import EarthRepeatOrbits

## following method to import all the functions from a module is not recommended 
## in the cheat sheets
#from Modules.BasicAstrodynamics import * 


TestCase = 2 # 1 - Coordinate conversions
             # 2 - Earth repeat orbits
             # 3 - Integ1 problem - Car motion - NI using Integrators module
             
############################
####### Coordinate conversions ###
############################

if TestCase == 1:

    # Import required packages    
    import math 
    import numpy as np 
    import spiceypy as spice
    import pygmo

    # Part(Problem) 1 of Basics-I assignment 
    S_bar = np.array([8751268.4691, -7041314.6869, 4846546.9938, 332.2601039, -2977.0815768, -4869.8462227]) 
    
    
    # Extracting mu from SPICE kernels
    print (spice.tkvrsn('TOOLKIT'))
    
    spice.furnsh("./External_files/Spice_kernels/kernel_load.txt")
    
    muE = spice.bodvrd( 'EARTH', 'GM', 1 )
    
    #mu = 398600.441E9 # Gravitational parameter for Earth [m^3/s^2]
    mu = muE[1][0] * 1e9
    
    spice.kclear
    
    CovertedKep = convertCartesianToKepler(S_bar,mu,True,True)    # Position arguments are passed
    
    # Part(Problem) 2 of Basics-I assignment 
    a = 12158817.9615 # Semi-major axis[m]
    e = 0.014074320051 # Eccentricity
    i = 52.666016957 # Inclination [deg]
    RAAN = 323.089150643 # Right ascention of ascending node [deg]
    omega = 148.382589129 # Argument of pericenter [deg]
    M = 112.192638384 # Mean anomaly[deg] 
    
    Kepler = np.array([a,e,i,RAAN,omega,M])
    
    ConvertedCarte = convertKeplerToCartesian(Kepler,mu,1,isInputInDegree = True, isPrint=True)
    ConvertedCarte2 = convertKeplerToCartesian(Kepler,mu,7,isInputInDegree = True, isPrint=True)
    
elif TestCase == 2:
    import numpy as np2
    
    Run = 3 # 1 - Approach 1 with i as unknown
            # 2 - Approach 2 with i as unknown 
            # 3 - Approach 2 with altitude as unknown

    # Approach 1, i as unknown 
    if Run == 1:
        jk = np2.matrix('39, 3; 40, 3; 41, 3; 42, 3; 43,3; 44, 3 ; 45,3; 46, 3; 47,3 ; 48, 3')
        e = 0
        a = np2.array([200,1200,1])
        Result = EarthRepeatOrbits(jk,e,a,'Alti',False,True)   
        
    # Approach 2, i as unknown 
    elif Run == 2:
        jk = np2.matrix('39, 3; 40, 3; 41, 3; 42, 3; 43,3; 44, 3 ; 45,3; 46, 3; 47,3 ; 48, 3')
        e = 0
        a = np2.array([200,1200,1])
        Result = EarthRepeatOrbits(jk,e,a,'Alti',True,True)   
        
    # Approach 2 with altitude as unknown
    elif Run == 3:
        jk = np2.matrix('14, 1; 43, 3; 29, 2; 59, 4; 74, 5 ; 15,1')
        e = 0
        i = np2.array([28,29,1])
        Result = EarthRepeatOrbits(jk,e,i,'Inclin',True,True)  
        
elif TestCase == 3:
    
    def CarProblem1D(S):
        "Defining the state derivative function for the car motion proble from Integ1"
        
        # error handling
        assert(S.size == 2),"Incorrect array size for the state variable passed to the STF."
        
    