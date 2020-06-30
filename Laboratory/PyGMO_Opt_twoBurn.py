# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Laboratory

Date of creation: Tue Jun 30 10:24:25 2020

Version: Base 

Description:
        1. 
        2. 
        3. 
        4. 
        5. 
        6. 
        7. Output validation - TBD [UPDATE THIS]
"""
###########################################
#### Importing the required modules  ######   
###########################################

import pygmo as pg
from Modules.AstrodynamicsProblems import TransferOrbits
import numpy as np

###########################################
## Definition: PyGMO problem ##############   
###########################################

class pagmo_TwoburnTransfer:
    
    # Creating a constructor, that creates an object of TransferOrbits class and set the proble parameters
    def __init__(self,para):
        # creating an instance of the "TransferOrbits" class
        self.__obj2BT = TransferOrbits(1)
        
        # Setting up parameter values
        self.parameters = para
        self.__obj2BT.setProblemParameters(para)
    
    # Defining the mandatory fitness function
    def fitness(self,x):
        DeltaVs = self.__obj2BT.computeDeltaV(x)
        return [sum(DeltaVs)]
         
    # Defining the mandatory get_bounds function    
    def get_bounds(self):
        return ([0],[self.parameters[1]])
   
#####################################
########### User inputs #############    
#####################################

# Parameters definition
para = np.array([185,28.5])

# Optimization parameters
nGen = 20


visualize = 2 # 0 - no visualizations
              # 1 - [UPDATE THIS]
              
######################################
######## Computations ################
######################################
problem = pg.problem(pagmo_TwoburnTransfer(para))
#print(problem)              

algo = pg.algorithm(pg.bee_colony(nGen, limit = 20))
pop = pg.population(problem,10)
pop = algo.evolve(pop)
print(pop.champion_f)               


#######################################
######### Visualization ###############
#######################################
#
#if visualize == 0:
#    print('No visualizations are requested by the user.')
#elif visualize == 1:
#    # [EDIT THIS]
#elif visualize ==2:
#    # [EDIT THIS]  
#else:
#    print('Unrecognized input for the visualization option.')

