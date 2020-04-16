# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Laboratory

Date of creation: Thu Apr 16 09:35:16 2020

Version: Base

Description:
    1. Integ-1 assignment to solve problem of motion of car using numerical integrators
    2. 

"""

import numpy as npMain

def Analytical(step,t0,te,X0,V0,a):
    # importing the required modules
    import numpy as np
    
    # Error handling
    assert ((te-t0)%step ==0), "Integer number of steps are not possible. Check the inputs!"
    
    # Creating variables of interest    
    time = np.arange(t0, te+step, step) # time vector
    AnaSol = np.zeros(time.size*3)
    AnaSol.shape = [time.size, 3] # Storage for the solution
    
    # Solving analytical expression at each time-step
    for count in range (time.size):
        AnaSol[count,0] = time[count]
        t = AnaSol[count,0]
        AnaSol[count,1] = X0*t + 0.5 *a* t**2 
        AnaSol[count,2] = V0 + a*t 
         
    return AnaSol 

def Euler(step,t0,te,S,a):
    # importing the required modules
    import numpy as np
    
    # Creating variables of interest
    time  = np.arange(t0,te+step,step) # time vector
    EulerSol = np.zeros(time.size*3) 
    EulerSol.shape = [time.size,3]
    EulerSol[0,0:3] = [time[0],S[0],S[1]] # initial state

    # Using Euler integration formulation at each step
    for count in range (1, time.size):
         S = EulerSol[count-1,1:3] # creates another array representing state at that epoch
#         S_dot = [S[1],a] # would lead to creation of list
         S_dot2 = np.array([S[1],a]) # correct way - derivative vector
         S = S + S_dot2 * step
         EulerSol[count,0:3] = [time[count],S[0],S[1]] # state at that epoch
    
    return EulerSol


def RK4(step,t0,te,S0,a):
    # importing the required modules
    import numpy as np
    
    # Creating variables of interest
    time = np.arange(t0,te+step,step) # time vector
    RK4Sol = np.zeros(time.size*3) 
    RK4Sol.shape = [time.size,3] # storage for solution
    RK4Sol[0,0:3] = [time[0],S0[0],S0[1]] #initialization
    
    # Using RK4 integration formulation at each step
    for count in range (1,time.size):
        t1 = time[count-1] # last epoch
        S0 = RK4Sol[count-1,1:3]
        K1 = np.array([S0[1],a]) 
     
        t2 = t1+ step/2
        S02 = S0 + step* 0.5 * K1
        K2 = np.array([S02[1],a])
            
        t3 = t1 + step/2
        S03 = S0 + step* 0.5 * K2 
        K3 = np.array([S03[1],a])
     
        t4 = t1 + step
        S04 = S0 + step * K3 
        K4 = np.array([S04[1],a])
        
        phi = (K1 + 2*K2 + 2*K3 + K4)/6
        S = S0 + phi * step
        RK4Sol[count,0:3] = [time[count],S[0],S[1]]
    
    return RK4Sol
    

# Parameters definition
step = 4  # [s]
t0 = 0
te= 60
X0 = 0 # [m]
V0 = 0 # [m/s]
S0 = npMain.array([X0,V0]) # intial state
a = 2 # [m/s^2]

# Solving the problem using equations of motion - analytical way
#AnaSol = Analytical(step,t0,te,X0,V0,a)

# Solving the problem using Euler integrator
step = 0.1
EulerSol = Euler(step,t0,te,S0,a)

# Solving the problem using RK4 integrator
step = 0.1
RK4Sol = RK4(step,t0,te,S0,a)

