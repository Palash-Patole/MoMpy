# -*- coding: utf-8 -*-
"""
: Sat Mar 25 2020, ‏‎08:30:56

version: Version 1

Description:
    1. Developed a function to compute Earth repeat orbits geometry when
        a or i, j, k, e are provided as inputs
    2. Low fidelity approach (default) where Earth-repeat orbit solutions are 
        computed cosindering the effect of $J_2$ perturbation on $\Omega$ only  
    3. High fidelity approach where Earth-repeat orbit solutions are 
        computed cosindering the effect of $J_2$ perturbation on $\Omega$, 
        $\omega$, and $M$ 
    4. Either altitude or inclination values (minimum,maximum,step size) are 
        to be provided.
    5. Low fidelity approach does not work when the inclination is known 
        while altitude is unknown
    6. For eccentric orbits, altitude values (input and output) should 
        corrospond to the perigee altitude.        
        
"""

import numpy as np2

def EarthRepeatOrbits(jk,e,Variable,VarType,isHighFidelity=False,printStatus=False):
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
    spice.furnsh("./External_files/Spice_kernels/kernel_load.txt")
    muE = spice.bodvrd('Earth','GM',1)
    mu = muE[1][0] #km3/s2 for earth 
    J2 = 1082.63E-6 #J2 for earth - import from spice
    RE = spice.bodvrd('EARTH','RADII',3)
    Re= RE[1][0]
    De = 86164.1004 # s, Sidereal day, to be extracted from spice
    k2 = 0.75 * J2 * math.sqrt(mu)* Re*Re
    r2d = 1/spice.rpd() # Radian to degree conversion
    
    # Creating storage
    steps = int((Variable[1]-Variable[0])/Variable[2])+1
    Result = np.zeros(int(jkSize[0])*int(steps)*6)
    Result.shape= [int(jkSize[0]),int(steps),6]
    
    # Computations
    for count in range(0,jkSize[0]):
#        print(count)
        for rowCount in range(0,int(steps)):
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
            
            # computations
            if isHighFidelity == False:
#                print('Approach 1')
                if VarType == 'Alti':
#                    print('Inclination is unknown')
#                    if e == 0:
#                        a = Re + var
#                        T = 2 * math.pi * math.sqrt(a**3/mu)
#                    else:
#                        print('Problem not defined yet for non-circular orbits.')
#                        sys.exit # DOES NOT WORK PROPERLY, LOOP IS NOT BROKEN
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
                        i = math.acos(DeltaL2/C1) * r2d # for i = [90, 180] deg 
                        Result[count,rowCount,4] = i
                        Nsolutions += 1
                    elif ( C3 < 0) & (abs(C3/C1) <=1): # this was missing in MATLAB script
                        DeltaL2 = C3
                        i = math.acos(DeltaL2/C1) * r2d # for i = [0, 90] deg 
                        Result[count,rowCount,5] = i   
                        Nsolutions += 1
                    elif (C3>0) & (abs(C3/C1) <=1):
                        DeltaL2 = C3
                        i = math.acos(DeltaL2/C1) * r2d # for i = [90, 180] deg 
                        Result[count,rowCount,5] = i   
                        Nsolutions += 1
                elif VarType == 'Inclin':
                    print('Function is not defined for low fidelity + unknown inclination case.')
                    return Result
            else:
#                print('Approach 2')
                if VarType == 'Alti':
#                    if e == 0:
#                        a = Re + var
#                    else:
#                        print('Problem not defined yet for non-circular orbits.')
#                        sys.exit
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
                        #       print(RAAN_dot)
                        omega_dot = k2 * math.pow(a0,-3.5) * (5*(math.cos(i/r2d))**2 -1) * math.pow(1-e**2,-2) * De * r2d
                        #       print(omega_dot)
                        M_dot = k2 * math.pow(a0,-3.5) * (3*(math.cos(i/r2d))**2 -1) * math.pow(1-e**2,-1.5) * De * r2d
                        #       print(M_dot)
                        n = (j/k) * (L_dot - RAAN_dot) -(omega_dot + M_dot)
                        #       print(n)
                        a1 = (mu/(n/(De*r2d))**2)**(1/3)
                        #       print('a1 is:',a1)
           
#                      # print('At iteration:',iterations,' a1 was ',a1, 'a0 was', a0)
                        if iterations > 1000:
#                            print('Maximum number of iterations executed')
                            a0 = a1
                            break
                        elif (abs(a1-a0) < 10**(-10)):  
#                            print('Convergence reached after ',iterations, 'iterations')    
                            a0 = a1
                            break
                        else:
                            a0 = a1    
                    
                    # Only tracking the feasible solutions
                    if ((a0 * (1-e)) - Re) > 0:
                        Nsolutions += 1
                        Result[count,rowCount,4] = (a0 * (1-e)) - Re
                    
    # Printing the status of solutions obtained
    if printStatus == True:
        print(Nsolutions, ' solutions are obtained for the Earth-repeat orbits.')
                                
            
#    print(Result.shape)
    return Result
    
    
##########################
######### Executing runs
##########################

Run = 1 # 1 - Approach 1 with i as unknown
        # 2 - Approach 2 with i as unknown 
        # 3 - Approach 2 with altitude as unknown

# Approach 1, i as unknown 
if Run == 1:
#    jk = np.matrix('14, 1; 43, 3; 29, 2; 59, 4; 74, 5 ; 15,1')
    jk = np2.matrix('39, 3; 40, 3; 41, 3; 42, 3; 43,3; 44, 3 ; 45,3; 46, 3; 47,3 ; 48, 3')
    e = 0.1
    a = np2.array([200,2000,5])
    Result = EarthRepeatOrbits(jk,e,a,'Alti',False,True)    

# Approach 2, i as unknown 
elif Run == 2:
    jk = np2.matrix('39, 3; 40, 3; 41, 3; 42, 3; 43,3; 44, 3 ; 45,3; 46, 3; 47,3 ; 48, 3')
    e = 0.1
    a = np2.array([200,2000,5])
    Result = EarthRepeatOrbits(jk,e,a,'Alti',True,True)   
    
# Approach 2 with altitude as unknown
elif Run == 3:
    jk = np2.matrix('14, 1; 43, 3; 29, 2; 59, 4; 74, 5 ; 15,1')
    e = 0
    a = np2.array([28,29,1])
    Result = EarthRepeatOrbits(jk,e,a,'Inclin',True,True)   