# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Modules

Date of creation: Tue May  5 02:40:00 2020

Version: developed with version 3 of the laboratory script, Integ3.py

Description:
        1. Definition of the base or superclass for numerical integrator
        2. Definition of the subclass - Euler integrator
        3. Definition of the subclass - RK4 integrator

"""
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
        assert(self.te>=self.t0),"Final epoch of integration is smaller than intial epoch."
        
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
             
            t2 = t1+ self.tStep/2
            S2 = S + self.tStep* 0.5 * K1
            K2 = self._sdf(S2)
            
            t3 = t1 + self.tStep/2
            S3 = S + self.tStep* 0.5 * K2 
            K3 = self._sdf(S3)
     
            t4 = t1 + self.tStep
            S4 = S + self.tStep * K3 
            K4 = self._sdf(S4)
        
            phi = (K1 + 2*K2 + 2*K3 + K4)/6
            S = S + phi * self.tStep
            for index in range(1,self.stateSize):
                Sol[count,index] = S[index-1] # saving current state
        
        return Sol