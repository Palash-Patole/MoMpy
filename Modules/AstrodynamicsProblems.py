# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Modules

Date of creation: Wed Jun  3 16:09:06 2020

Version: Base 

Description:
        1. Module to contain useful astrodynamic problems
        2. Circular LEO to GEO two-burn transfer problem - Problem type 1
        3. Circular LEO to another circular-equatorial orbit at higher altitude - Problem type 2
        4. Two burn transfer problem verfied by checking the function and its derivative value plot. Comparing the known results for minima.
        5. 
        6. 
        7. 
"""
###############################################################
#### Importing the required modules - for local testing  ######   
###############################################################

import numpy as np
import matplotlib.pyplot as plt
import math


###########################################
## Definition: Two-burn trasfer problem ###   
###########################################
class TransferOrbits:
    """ Defining two/three burn transfer problem from/to LEO to/from GEO/Higher altitude orbit """
    
    # constructor - assigning and checking the problem type
    def __init__(self,problemtype):
        self.problemType = problemtype 
        self.__paraSet = False
        
        # error handling
        assert(self.problemType==1 or self.problemType==2),"Currently, only problem type 1 and 2 are supported. Check the descriptions of supported problems using 'stateProblemTypes' method of the class."

        # setting up the number of required parameters, to be provided by user
        if self.problemType==1:
            self.__requiredPara = 2
        elif self.problemType == 2:
            self.__requiredPara = 4
    
    # Setting up the parameters of the problem
    def setProblemParameters(self,ParaSet):   
        self.__Paraset = ParaSet
        
        # Error-handling
        assert(self.__Paraset.size==self.__requiredPara),"Incorrect array-size for the parameter matrix."
        
        
        # Extracting required parameters using SPICE kernels if possible
        try:
            import spiceypy as spice
            import sys

            sys.path.append('../') # add everything that exists in the MoMPy directory
            
            spice.furnsh('../External_files/Spice_kernels/Load_Kernels3.txt')
            muE = spice.bodvrd('EARTH','GM',1)
            mu = muE[1][0] # mu for Earth [km^3/s^2]
            RE = spice.bodvrd('EARTH','RADII',3)
            Re= RE[1][0] # [km], Average radius of Earth
            spice.kclear    
        except:
            mu = 398600.435436095925979316234588623046875 # manually entering in case the import fails
            Re = 6378.136599999999816645868122577667236328125 # manually entering in case the import fails
        finally:
            hGEO = 35768 # Altitude of GEO [km]
            
        # Importing required modules
        import math
        
        # Defining/setting parameters that are specific to the problem-type
        if self.problemType==1:
            # Extraction
            hLEO = self.__Paraset[0] # Altitude of circular LEO [km]
            self.__iLEO = self.__Paraset[1] # Inclination of LEO [deg]
            assert( hLEO < hGEO ), "LEO altitude can not be greater than GEO altitude. Check inputs."
            
            # Computation
            ra = Re + hGEO # apoapsis of the GTO [km]
            rp = Re + hLEO # periapsis of the GTO [km]
            self.__VLEO = math.sqrt(mu/rp)  # [km/s]
            self.__VGEO = math.sqrt(mu/ra) # [km/s]          
            e = (ra-rp)/(ra+rp) # eccentricity of the GTO
            a = (ra +rp)/2 # semi-major axis of the GTO
            self.__Va = math.sqrt(mu * (1-e)/(a * (1+e))) # Apoapsis velocity in the GTO [km/s] 
            self.__Vp = math.sqrt(mu * (1+e)/(a * (1-e))) # Periapsis velocity in the GTO [km/s]
        
            self.__paraSet = True # Parameters are now completely set and other methods can be used
        elif self.problemType == 2:
            # Extraction
            self.__r1LEO = Re + self.__Paraset[0] # Semi-major axis of the circular LEO [km]
            self.__iLEO = self.__Paraset[1] # Inclination of the circular LEO [deg]
            self.__r2HEO = Re + self.__Paraset[2] # Semi-major axis of the high altitude final orbit [km]
            self.__raInterMedLimit = Re + self.__Paraset[3] # Limiting value of the apoapsis distance of the intermediate transfer orbit [km]
            assert( self.__r2HEO < self.__raInterMedLimit ), "Limiting value of the apoapsis distance of the intermediate transfer orbit has to be greater than semi-major axis of the final circular-equatorial orbit. Check inputs."
            
            # Computation
            self.__VC1 = math.sqrt(mu/self.__r1LEO) # Velocity of the circular LEO [km/s]
            self.__VC2 = math.sqrt(mu/self.__r2HEO) # Velocity of the circular final orbit [km/s]            
            
            self.__paraSet = True # Parameters are now completely set and other methods can be used
        else:
            print('Problem type is not supported.')
        
    # Method to compute Delta Vs for the given set of user inputs
    def computeDeltaV(self,variablesSet):
        
        # Error handling
        assert(self.__paraSet == True),"First set the parameters of the object using the method - 'setProblemParameters'."
        
        # importing the required modules
        import math
        import numpy as np
        
        if self.problemType==1:
            # storage for the results
            DeltaVs = np.zeros(2)
            
            # computations
            i1 = variablesSet * math.pi/180
            i2 = (self.__iLEO - variablesSet) * math.pi/180
            
            DeltaVs[0] = math.sqrt(self.__VLEO**2 + self.__Vp**2 - 2*self.__VLEO*self.__Vp*math.cos(i1)) # DeltaV for first burn
            DeltaVs[1] = math.sqrt(self.__VGEO**2 + self.__Va**2 - 2*self.__VGEO*self.__Va*math.cos(i2)) # DeltaV for second burn
        
        elif self.problemType==2:
            # storage for the results
            DeltaVs = np.zeros(3)
            
            # computations
            i2 = variablesSet[0] * math.pi/180 # Inclination change at the second burn
            i3 = (self.__iLEO - variablesSet[0]) * math.pi/180 # Inclination change at the third burn
            ra = variablesSet[1] * (self.__raInterMedLimit-self.__r2HEO) # Apoapsis altitude of the intermediate transfer orbit [km]
                        
            et1 = (ra-self.__r1LEO)/(ra + self.__r1LEO) # Eccentricity of the first transfer leg 
            et2 = (ra-self.__r2HEO)/(ra + self.__r2HEO) # Eccentricity of the second transfer leg
            VCa = self.__VC1 * math.sqrt(self.__r1LEO/ra) # Velocity of the circular orbit with radius same as the apoapsis distance of the intermediate transfer orbit
            
            # First Delta-V required, to go from circular LEO to first leg of the intermediate elliptical transfer orbit, orbital plane is not changed 
            DeltaVs[0] = self.__VC1 * (math.sqrt(1+et1)-1)
            # Second Delta-V required, at the apoapsis of the first leg of the intermediate ellitptical transfer orbit, to put S/C on the second leg, an inclination change is achieved
            DeltaVs[1] = VCa * math.sqrt( (1-et1) + (1-et2) -2 * math.sqrt((1-et1)*(1-et2)) * math.cos(i2) )
            # Third Delta-V required, at the periapsis of the second leg of the intermediate ellitptical transfer orbit, to put S/C in the final circular orbit, an inclination change is achieved
#            DeltaVs[2] = self.__VC2 * (math.sqrt(1+et2)-1) #Original formula in wakker when there is no inclination change at this maneuver
            DeltaVs[2] = self.__VC2 * math.sqrt( (1+et2) + 1 - 2 *math.sqrt(1+et2) * math.cos(i3))
        else:
            print('Problem type not supported.')
            
        return DeltaVs
    
    # Method to compute the required derivative value
    def computeDerivative(self,variablesSet):    
        # Error handling
        assert(self.__paraSet == True),"First set the parameters of the object using the method - 'setProblemParameters'."
        
        # importing the required modules
        import math
        
        if self.problemType==1:
            # computations
            i1 = variablesSet * math.pi/180
            i2 = (self.__iLEO - variablesSet) * math.pi/180
            
            A = self.__VLEO * self.__Vp * math.sin(i1) / math.sqrt(self.__VLEO**2 + self.__Vp**2 - 2*self.__VLEO*self.__Vp*math.cos(i1))
            B = self.__VGEO * self.__Va * math.sin(i2) / math.sqrt(self.__VGEO**2 + self.__Va**2 - 2*self.__VGEO*self.__Va*math.cos(i2))
            
            return (A-B)
    
    # Documenting problem types
    @staticmethod 
    def stateProblemTypes():
        problems = ["Problem type 1: Circular LEO to GEO, two-burn transfer", "Problem type 2: Circular LEO to Higher altitude circular-equatorial orbit, three-burn transfer"]
        return problems

              
######################################
######## Local testing ###############
######################################
        
problem = 1 # Problem to be solved
            # 0 - No local testing
            # 1 - two-burn transfer problem, check inputs below
            # 2 - three burn transfer problem 
visualize = False # When locally tested, visualizes the validation results            

if problem==1:
# Setting up the two-burn transfer problem
    object1 = TransferOrbits(1)
    para = np.array([185,28.5])
    object1.setProblemParameters(para)
    step = 0.001
    
    # Storage
    i1Data = np.arange(0,para[1]+0.1,step)
    Results = np.zeros(i1Data.size*5)
    Results.shape = [i1Data.size,5]
    
    # Solving problem at various i1 values
    rowIndex = 0
    for i1 in i1Data:
        DeltaVs = object1.computeDeltaV(i1)
        Results[rowIndex,0] = i1
        Results[rowIndex,1] = DeltaVs[0]+DeltaVs[1]
        Results[rowIndex,2] = DeltaVs[0]
        Results[rowIndex,3] = DeltaVs[1]
        Results[rowIndex,4] = object1.computeDerivative(i1)
        rowIndex +=1
    
    # finding row index where the functions's value is minimum    
    FminColumn = np.where(Results[:,1]==np.amin(Results[:,1])) 
    FminColumn = FminColumn[0][0]
    # finding index value where function's derivative has a zero value
    dFminColumn = np.where(abs(Results[:,4])==np.amin(abs(Results[:,4])))
    dFminColumn == dFminColumn[0][0]
    
    # Checking where is minima is found and where it is expected based on the derivative function value:
    if (FminColumn==dFminColumn):
        print("The function has minimum value where the derivative function has zero/closest to zero value.")
    
    print('\n\nThe minimum function is '+ str(Results[FminColumn,1]) + ' and occurs at '+ '\u0394i1 :'+ str(Results[FminColumn,0]) +'.')
    print('The derivative function has a value closest to zero: '+str(Results[dFminColumn,4])+' and it occurs at '+ '\u0394i1 :'+ str(Results[dFminColumn,0]) +'.')
    
elif problem==2:    
    # For generating first plot of Figure 13.9 (pp324), from Wakker Astrodynamics notes
    object2 = TransferOrbits(2)
    
    mvalues = np.array([2,5,10,30]) # m-values (ra_transfer/rLEO) to be tested
    altitudeLEO = 185 # Altitude of the intial LEO orbit [km]
    inclinLEO = 50 # Inclination of the LEO orbit [deg]
    muE = 398600.435436095925979316234588623046875 # GM for Earth [km3/s2]
    RE = 6378.136599999999816645868122577667236328125 # Average radius of Earth [km]
    limitingDistance = 1000000 # Maximum apoapsis distance of the intermediate transfer orbit [km]
    Vc1 = math.sqrt(muE/(RE+altitudeLEO)) # Velocity in LEO orbit [km/s], to be used for normalization of results
    
    # Creating empty storage for results various values of m
    mResults = np.zeros(3)
    # Tracking the number of solutions obtained
    solTracker = np.zeros(mvalues.size*2)
    solTracker.shape = [mvalues.size,2]
    # Computing solutions based on the formulas in Wakker
    DVfromformula = np.zeros(6)
    # computations and various n and m values
    for m in mvalues:
        solTracker[np.where(mvalues==m),0] = m
        for n in np.arange(0,21.0,0.1): # n-values (rHEO/rLEO) to be tested
            if (m > n) & (n >0):
                # Knowing the altitudes/distances of interest
                hHEO = n*(RE+altitudeLEO)-RE # Altitude of the final orbit [km]
                ra = m * (RE+altitudeLEO) # Apoapsis distance of the intermediate transfer orbit [km]
                
                
                # Setting up the three-burn transfer problem
                para2 = np.array([altitudeLEO,inclinLEO,hHEO,limitingDistance-RE])
                object2.setProblemParameters(para2)
                
                # Computing for the rHEO distance
                fraction =  ra / (limitingDistance-(hHEO+RE))
                variables = np.array([inclinLEO,fraction])
                DeltaVs2 = object2.computeDeltaV(variables)

                # Storing the result
                normaltotalDeltaV = np.sum(DeltaVs2)/Vc1 
                mResults = np.vstack((mResults,np.array([n,m,normaltotalDeltaV])))
                
                # Computing the Delta V based on formulas in Wakker
                nDV1 = math.sqrt(2*m/(m+1)) - 1 # DeltaV1/Vc1
                nDV2 = math.sqrt((2/m) * (1/(m+1) + n/(m+n) - 2 * math.sqrt(n/((m+1) * (m+n)))* math.cos(inclinLEO*math.pi/180)  ) ) # DeltaV2/Vc1
                nDV3 = math.sqrt(1/n) * (math.sqrt(2*m/(m+n)) - 1) # DeltaV3/Vc1
                
                # Storing these formula-based results
                DVfromformula = np.vstack((DVfromformula,np.array([n,m,nDV1,nDV2,nDV3,nDV1+nDV2+nDV3])))
                
                
        solTracker[np.where(mvalues==m),1] = (mResults.size/3)-1 # tracking how many solutions were found for this m value, to visualize it later
    
    mResults = np.delete(mResults,0,0) # Deleting the first empty row
    DVfromformula = np.delete(DVfromformula,0,0) # Deleting the first empty row
    
######################################
######## Visualization ###############
######################################   
    


if (visualize==True):    
    if problem==1:          
        # Double axis plot of function and its derivative value
        fig, ax1 = plt.subplots()
        
        # First axis
        color  = 'tab:blue'
        ax1.set_xlabel('First-maneuver inclination change '+r'$(\Delta i_1)$'+' [deg]')
        ax1.set_ylabel('Total '+r'$\Delta V$' + ' [m/s]')
        ax1.plot(Results[:,0],Results[:,1],color=color)
        ax1.tick_params(axis='y',labelcolor=color)
        ax1.set_xlim(left=i1Data[0])
        ax1.set_xlim(right=i1Data[-1])
        plt.axvline(x=Results[FminColumn,0],color='k',linestyle='--',label="Function has minimum value here")
        plt.axvline(x=Results[dFminColumn,0],color='g',linestyle='-.',label="Derivative has zero value here") 
        plt.legend()
            
        ax2 = ax1.twinx() # second y-axis
        color = 'tab:red'
        ax2.set_ylabel('Derivative of total '+ r'$\Delta V$'+ ' function')
        ax2.plot(Results[:,0],Results[:,4],color=color)
        
        ax2.tick_params(axis='y',labelcolor=color)
        
        #common properties for the figure
        fig.tight_layout()
        plt.title('Two-burn trasfer problem: function and its derivative value')
        plt.grid(True)
        plt.show()
        
        
    if problem==2:
        fig = plt.figure()
        fig.suptitle("Validation of results obtained with 'TransferOrbits' class for three-burn transfer problem", fontsize=16)
        plt.subplot(221)
        plt.plot(mResults[0:int(solTracker[0,1]),0],mResults[0:int(solTracker[0,1]),2],'o-',label='Using MoMPy class, m=2')
        plt.plot(mResults[0:int(solTracker[0,1]),0],DVfromformula[0:int(solTracker[0,1]),5],'*-',label='Using formulae from Wakker,m=2')
        plt.xlim(left = 0)
        plt.xlabel('n')
        plt.ylabel('Total '+ r'$\Delta V/ V_{c_1}$')
        plt.grid(True)    
        plt.legend()
        
        plt.subplot(222)
        plt.plot(mResults[int(solTracker[0,1]):int(solTracker[1,1]),0],mResults[int(solTracker[0,1]):int(solTracker[1,1]),2],'o-',label='Using MoMPy class, m=5')    
        plt.plot(mResults[int(solTracker[0,1]):int(solTracker[1,1]),0],DVfromformula[int(solTracker[0,1]):int(solTracker[1,1]),5],'*-',label='Using formulae from Wakker,m=5')
        plt.xlim(left = 0)
        plt.xlabel('n')
        plt.ylabel('Total '+ r'$\Delta V/ V_{c_1}$')
        plt.grid(True)    
        plt.legend()
        
        plt.subplot(223)
        plt.plot(mResults[int(solTracker[1,1]):int(solTracker[2,1]),0],mResults[int(solTracker[1,1]):int(solTracker[2,1]),2],'o-',label='Using MoMPy class, m=10')
        plt.plot(mResults[int(solTracker[1,1]):int(solTracker[2,1]),0],DVfromformula[int(solTracker[1,1]):int(solTracker[2,1]),5],'*-',label='Using formulae from Wakker,m=10')
        plt.xlim(left = 0)
        plt.xlabel('n')
        plt.ylabel('Total '+ r'$\Delta V/ V_{c_1}$')
        plt.grid(True)    
        plt.legend()
        
        plt.subplot(224)
        plt.plot(mResults[int(solTracker[2,1]):int(solTracker[3,1]),0],mResults[int(solTracker[2,1]):int(solTracker[3,1]),2],'o-',label='Using MoMPy class, m=30')
        plt.plot(mResults[int(solTracker[2,1]):int(solTracker[3,1]),0],DVfromformula[int(solTracker[2,1]):int(solTracker[3,1]),5],'*-',label='Using formulae from Wakker,m=30')
        plt.xlim(left = 0)
        plt.xlabel('n')
        plt.ylabel('Total '+ r'$\Delta V/ V_{c_1}$')
        plt.grid(True)    
        plt.legend()