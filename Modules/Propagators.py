# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Modules

Date of creation: Thu May  7 06:47:10 2020

Version: developed with version 3 of the laboratory script, Integ3.py

Description:
        1.Defined the propagator class 
        2.Propagators supported so far: Cowell
        3.Central bodies supported so far: Earth
        4.Definted the pointMassGravity function to compute state derivative with Cowell propagator  
"""
###########################################
### Definition: Propagator Class ##########   
###########################################
    

class Propagator():
    """ Definition of the propagator class""" 
    
    # constructor
    def __init__(self,propagator):
        self.propagator = propagator
        self.centralBodyName = ""
        self.__propagatorSet = False
        self.__centralBodySet = False
        self.__stateSizerequired = 0
        
        # error handling
        if self.propagator == "Cowell": 
            self.__stateSizerequired = 6
            self.__propagatorSet = True
        else:
            import warnings
            warnings.warn('Specified propagator is not supported yet or an incorrect propagator name. Use listAllpropagators() class method to know the supported propagator names.')
            
    # documenting the supported propagators type            
    @staticmethod 
    def listAllpropagators():
        supportedPropagators = ["Cowell"]
        return supportedPropagators
    
    # setting up the central body         
    def setCentralBody(self,CentralBN):
        self.centralBodyName = CentralBN 
        
        ## fix following to extract the mu from SPICE
        #{
        # import the required modules
#        import spiceypy as spice
#        import sys
        
#        sys.path.append('../') # add everything that exists in MoMpy directory
#        
#        # extract the gravitional parameter from the SPICE kernels + error handling
#        spice.furnsh("../External_files/Spice_kernels/Load_Kernels2.txt")
    
        #}
        
        if (self.centralBodyName == "Earth"):
#           # FIX following when SPICE kernels are loaded sucessfully 
#            CBmu = spice.bodvrd('EARTH','GM',1)
#            self.__CBmu = CBmu[1][0]*1e9 # Gravitational parameter for Earth [m^3/s^2]
#            spice.kclear
            
            self.__CBmu = 398600435436095.94 # mu for Earth
            self.__centralBodySet = True
        else:
            import warnings
            warnings.warn('Name of the central body is not supported yet. Use listAllCentralBodies() class method to know all the celestial bodies supported as a central body.')    
            self.__CBmu = None
    
    # documenting the central bodies supported
    @staticmethod
    def listAllCentralBodies():
        supportedCentralbodies = ["Earth"]
        return supportedCentralbodies

    ###########################################
    # Definition: Point mass gravity function #   
    ###########################################

    def pointMassGravity(self,S):
        "Function to compute point mass gravity acceleration and derivate function value for a state of the S/C"
        
        # error handling
        assert(self.__propagatorSet == True),"Invalid propagator selection."
        assert(self.__centralBodySet == True),"Invalid specification for the central body."
        assert(S.size == self.__stateSizerequired),"Incorrect size for the state argument of the pointMassGravity function"
    
        # importing required modules
        import math
        import numpy as np
    
        # computing the acceleration
        if (self.propagator == 'Cowell'):
            r = math.sqrt(S[0]**2 + S[1]**2 + S[2]**2)
            S_dot = np.array([S[3], S[4], S[5], -self.__CBmu*S[0]/(r**3), -self.__CBmu*S[1]/(r**3), -self.__CBmu*S[2]/(r**3)])
        
        return S_dot
            
###########################################
## Testing locally - class and object #####  
###########################################            

#object1 = Propagator("Cowell")
#object1.setCentralBody("Earth")
#print(Propagator.listAllCentralBodies())
#print(Propagator.listAllpropagators())



