# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: (base folder)
Date of creation: Fri Nov  8 20:57:12 2019

version: Base

Description:
    1. Generic canvas to test different modules
    2. Tested - coordinate coversations
    3. 3 options of computing the Earth-repeat orbit solutions
    4. 
"""

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
    from Modules.BasicAstrodynamics import convertCartesianToKepler
    from Modules.BasicAstrodynamics import convertKeplerToCartesian

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
    # import the required modules
    import numpy as np2
    from Modules.OrbitDesign import EarthRepeatOrbits
    import pickle # to export the variable data
    
    Run = 3 # 1 - Approach 1 with i as unknown
            # 2 - Approach 2 with i as unknown 
            # 3 - Approach 2 with altitude as unknown

    # Approach 1, i as unknown 
    if Run == 1:
        jk = np2.matrix('39, 3; 40, 3; 41, 3; 42, 3; 43,3; 44, 3 ; 45,3; 46, 3; 47,3 ; 48, 3')
        e = 0
        a = np2.array([200,1200,1])
        Result = EarthRepeatOrbits(jk,e,a,'Alti',False,True) 
        
        # Saving the results for post-processing
        f = open('Output/EarthRepeats_data1.pckl', 'wb')
        pickle.dump(Result, f)
        f.close()
        
    # Approach 2, i as unknown 
    elif Run == 2:
        jk = np2.matrix('39, 3; 40, 3; 41, 3; 42, 3; 43,3; 44, 3 ; 45,3; 46, 3; 47,3 ; 48, 3')
        e = 0
        a = np2.array([200,1200,1])
        Result = EarthRepeatOrbits(jk,e,a,'Alti',True,True) 
        
        # Saving the results for post-processing
        f = open('Output/EarthRepeats_data2.pckl', 'wb')
        pickle.dump(Result, f)
        f.close()
        
    # Approach 2 with altitude as unknown
    elif Run == 3:
        jk = np2.matrix('14, 1; 43, 3; 29, 2; 59, 4; 74, 5 ; 15,1')
        e = 0
        i = np2.array([28,29,1])
        Result = EarthRepeatOrbits(jk,e,i,'Inclin',True,True)  
        
        # Saving the results for post-processing
        f = open('Output/EarthRepeats_data3.pckl', 'wb')
        pickle.dump(Result, f)
        f.close()
        
elif TestCase == 3:
    
    # importing the required modules 
    import numpy as np    
    from Modules.Integrators import *
    from decimal import *
    import matplotlib.pyplot as plt
    
    # defining the state-derivative function
    def CarProblem1D(S):
        "Defining the state derivative function for the car motion proble from Integ1"
        
        # error handling
        assert(S.size == 2),"Incorrect array size for the state variable passed to the STF."
        
        # defining a constant acceleration
        a = 2 # [m/s^2]    
        
        Sdot = np.array([S[1],a])
        
        return Sdot
    
    # Parameters definition
    t0 = 0
    te= 60
    X0 = 0 # [m]
    V0 = 0 # [m/s]
    S0 = np.array([X0,V0]) # intial state
    steps = np.array([10, 1, 0.5, 0.1,0.01,0.001])
    ReferenceSol = 3600 # Reference (true) solution for the final distance, to compute % error    
    printRK4results = False # Prints the final distance computed using the RK4 integrator using a high precision output
    
    
    # Creating storage for the results of the step-size-variation analysis
    EulerResults = np.zeros(steps.size * 3) # saving step size, the number of function evaluations, and the final distance
    EulerResults.shape = [steps.size,3]
    
    # Incorrect way
#    RK4Results = EulerResults # same empty storage variable for both integrators 
    
    RK4Results = np.zeros(steps.size * 3) # saving step size, number of function evaluations, and the final distance
    RK4Results.shape = [steps.size,3]
    
    # Creating objects for the Euler and RK4 Integrator class
    object1 = EulerIntegrator(S0)
    object2 = RK4Integrator(S0)
    
    # Setting up state derivate function
    object1.setStateDerivativeFunction(CarProblem1D)
    object2.setStateDerivativeFunction(CarProblem1D)
    
    # Solving for different step-sizes 
    rowIndex  = 0
    for step in steps:
        time = np.array([t0,te,step])
        result1 = object1.integrate(time)
        result2 = object2.integrate(time)
        
        # Saving the step-size, the number of function evaluations, and final distance
        EulerResults[rowIndex,0] = step
        EulerResults[rowIndex,1] = len(result1) - 1
        EulerResults[rowIndex,2] = result1[-1,1]
        
        RK4Results[rowIndex,0] = step
        RK4Results[rowIndex,1] = len(result2) - 1
        RK4Results[rowIndex,2] = result2[-1,1]
        
        # Printing results with high precision, for the RK4 integrator
        if printRK4results == True :
            getcontext().prec = 24
            print('At step size of ',step, ' , final distance (RK4 integrator): ',Decimal(RK4Results[rowIndex,2]))
        
        rowIndex +=1 
        
    plt.semilogx(EulerResults[:,1],(ReferenceSol-EulerResults[:,2])*100/ReferenceSol,'*-',label='Euler Integrator')
    plt.semilogx(RK4Results[:,1],(ReferenceSol-RK4Results[:,2])*100/ReferenceSol,'+-',label='RK4 Integrator')
    plt.xticks(RK4Results[:,1],[str(RK4Results[0,1]),str(RK4Results[1,1]),str(RK4Results[2,1]),str(RK4Results[3,1]),str(RK4Results[4,1]),str(RK4Results[5,1])])
    plt.xlabel('Number of function evaluations')
    plt.ylabel('% error in the final computed distance')  
    plt.title('Accuracy of numerical integrators vs. number of function evaluations')
    plt.grid(True)
    plt.legend()
    
    