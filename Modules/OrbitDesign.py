# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Modules
Date of creation: Fri Mar 27 2020, 07:46:28 

version: Base

Culmination of:
    1. Ver2 of Basics-VII_ver1.py - 
       functions for finding Earth-repeat orbit solutions
       Added on: March 27, 2020
"""

def EarthRepeatOrbits(jk,e,Variable,VarType,isHighFidelity=False,printStatus=False):
    """Find a set of repeating ground track orbits around Earth (Earth-repeat orbits)."""
   
    # importing the required modules
    import spiceypy as spice
    import math
    import numpy as np
    
    Result = 0
    Nsolutions = 0
      
    # Error handling
    jkSize = jk.shape
    if (int(jkSize[1])!=2):
        print('Incorrect matrix size for jk matrix. Only two coloumns are required.')
        return Result
    elif (VarType != 'Alti') and (VarType != 'Inclin'):   
        print('Inrecognized input for the argument specifying the variable type.')
        return Result
    elif e >= 1.0 or e < 0:
        print('Only circular or elliptical orbits are possible. Check the value of eccentricity.')
        return Result
    elif (Variable[1]-Variable[0])% Variable[2] !=0:
        print('Integer number of steps are not possible. Check inputs for the argument - "Variable".')
        return Result
        
    # Extracting the parameters
    try:
        spice.furnsh("./External_files/Spice_kernels/Load_Kernels.txt")
        muE = spice.bodvrd('Earth','GM',1)
        mu = muE[1][0] # [km3/s2] for earth     
        RE = spice.bodvrd('EARTH','RADII',3)
        Re= RE[1][0] # [km], Average radius of Earth
        spice.kclear
    except:
        mu = 398600.435436095925979316234588623046875 # manually entering in case the import fails
        Re = 6378.136599999999816645868122577667236328125 # manually entering in case the import fails
        
    J2 = 1082.63E-6 #J2 for earth    
    De = 86164.1004 # [s], Sidereal day
    k2 = 0.75 * J2 * math.sqrt(mu)* Re*Re
    r2d = 1/spice.rpd() # Radian to degree conversion
    
    # Creating storage
    steps = int((Variable[1]-Variable[0])/Variable[2])+1
    Result = np.zeros(int(jkSize[0])*int(steps)*6)
    Result.shape= [int(jkSize[0]),int(steps),6]
    
    # Computations
    for count in range(0,jkSize[0]): #looping over j and k values
        
        for rowCount in range(0,int(steps)): #looping over variable values
            j = jk[count,0]    
            k = jk[count,1]    
            
            # Extracting the value of variable
            var = Variable[0]+rowCount*Variable[2]
            
            # Storing values known so far
            Result[count,rowCount,0] = j  
            Result[count,rowCount,1] = k  
            Result[count,rowCount,2] = e  
            Result[count,rowCount,3] = var  
            Result[count,rowCount,4] = math.nan # to be computed 
            Result[count,rowCount,5] = math.nan # to be computed
            
            if isHighFidelity == False:
                if VarType == 'Alti':
                    a = (Re + var)/(1-e)
                    T = 2 * math.pi * math.sqrt(a**3/mu)
                    
                    DeltaL1 = -2 * math.pi * (T/De)          
                    # DeltaL2 = C1 * cosi
                    C1 = (-3*math.pi * J2 * Re**2)/(a**2 * (1 - e**2)**2)
                    C2 = -2 * math.pi* k/j  - DeltaL1
                    C3 = 2 * math.pi* k/j  - DeltaL1
                    
                    # Inverse cosine
                    if (C2 <0) & (abs(C2/C1) <=1):
                        DeltaL2 = C2
                        i = math.acos(DeltaL2/C1) * r2d # for i = [0, 90] deg 
                        Result[count,rowCount,4] = i
                        Nsolutions += 1
                    elif (C2 > 0) & (abs(C2/C1) <=1):
                        DeltaL2 = C2
                        i = math.acos(DeltaL2/C1) * r2d # for i = (90, 180] deg 
                        Result[count,rowCount,4] = i
                        Nsolutions += 1
                    elif ( C3 < 0) & (abs(C3/C1) <=1): 
                        DeltaL2 = C3
                        i = math.acos(DeltaL2/C1) * r2d # for i = [0, 90] deg 
                        Result[count,rowCount,5] = i   
                        Nsolutions += 1
                    elif (C3>0) & (abs(C3/C1) <=1):
                        DeltaL2 = C3
                        i = math.acos(DeltaL2/C1) * r2d # for i = (90, 180] deg 
                        Result[count,rowCount,5] = i   
                        Nsolutions += 1
                elif VarType == 'Inclin':
                    print('Function is not defined for low fidelity + unknown inclination case.')
                    return Result
            else:
                if VarType == 'Alti':
                    a = (Re + var)/(1-e)    
                    C1 = -2 * k2 * (a**-3.5) * ((1-e**2)**-2) * De*r2d  # RAAN_dot = C1 * cos(i)
                    C2 = k2 * (a**-3.5) *((1-e**2)**-2) * De*r2d # omega_dot = C2 *(5*(cosd(i))^2 -1)
                    C3 = k2 * (a**-3.5) *((1-e**2)**-1.5) * De * r2d # M_dot = C3 *(3*(cosd(i))^2 -1)
                    n = De * r2d * math.sqrt(mu/a**3);
                    
                    # For the quadratic equation A*x^2 + B*x +C = 0 with x = cos(i) 
                    A = 5* C2 + 3 * C3
                    B = C1 * j/k
                    C = n - 360 *j/k - (C2 + C3)
                    D = B**2-4*A*C    
                    
                    if D >=0:
                        if (D == 0) & (abs(-B/(2*A)) <=1):
                            i1 = math.acos(-B/(2*A))
                            Result[count,rowCount,4] = i1*r2d
                            Nsolutions += 1
                        elif (D > 0): 
                            x1 = (-B + math.sqrt(D))/ (2 * A)
                            x2 = (-B - math.sqrt(D))/ (2 * A)
                            if abs(x1) <=1:
                                i1 = math.acos(x1)
                                Result[count,rowCount,4] = i1*r2d    
                                Nsolutions += 1
                            if abs(x2) <=1:
                                i2 = math.acos(x2) 
                                Result[count,rowCount,5] = i2*r2d 
                                Nsolutions += 1
                elif VarType == 'Inclin':

                    a0 = math.pow((mu*De**2*k**2/(4*math.pi**2*j**2)),1/3)
                    iterations = 0
                    i = var
                    L_dot = 360
       
                    while True:
                        iterations = iterations+1
           
                        RAAN_dot = - 2 * k2 * math.pow(a0,-3.5) * math.cos(i/r2d) * math.pow(1-e**2,-2) * De * r2d
                        omega_dot = k2 * math.pow(a0,-3.5) * (5*(math.cos(i/r2d))**2 -1) * math.pow(1-e**2,-2) * De * r2d
                        M_dot = k2 * math.pow(a0,-3.5) * (3*(math.cos(i/r2d))**2 -1) * math.pow(1-e**2,-1.5) * De * r2d
                        n = (j/k) * (L_dot - RAAN_dot) -(omega_dot + M_dot)
                        
                        a1 = (mu/(n/(De*r2d))**2)**(1/3)

                        if iterations > 1000:
                            a0 = a1
                            break
                        elif (abs(a1-a0) < 10**(-10)):  
                            a0 = a1
                            break
                        else:
                            a0 = a1    
                    
                    # Only tracking the feasible solutions i.e. positive altitudes
                    if ((a0 * (1-e)) - Re) > 0:
                        Nsolutions += 1
                        Result[count,rowCount,4] = (a0 * (1-e)) - Re

    # Printing the status of solutions obtained
    if printStatus == True:
        print(Nsolutions, ' solutions are obtained for the Earth-repeat orbits.')

    return Result