# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Laboratory

Date of creation: Tue Jun 30 10:24:25 2020

Version: Base 

Description:
        1. Optimization of the two-burn transfer problem from AstrodynamicsProblems/TransferOrbits
        2. Validation against the optimum solution found using grid search in AstrodynamicsProblems_ver0.py, problem = 1 [Local testing]
        3. 
        4. 
        5. 
        6. 
        7. Output validation - TBD [UPDATE THIS]
"""
###########################################
#### Importing the required modules  ######   
###########################################
import sys
sys.path.append('../')
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
refSol = 4.272327966969703
refX = 2.165

# Optimization parameters
nGen = 20

# Alorithms of interest
Algos = {1:'Differential evolution',2:'Artifical Bee Colony'}

visualize = False # TBC
              
######################################
######## Computations ################
######################################
# Creating a storage variable
Results = np.zeros((len(Algos),6))              

# Setting the problem
problem = pg.problem(pagmo_TwoburnTransfer(para))
             
for i in range(1,len(Algos)+1):
    print('\n\n\n')
    print("******************************************************************")
    print('The selected algorithm is: ', Algos[i])
    print("******************************************************************")
    if i==1:
        algo = pg.algorithm(pg.de(nGen))
    elif i==2:
        algo = pg.algorithm(pg.bee_colony(nGen, limit = 20))
    elif i==3:
        pass
    
    pop = pg.population(problem,10)
    pop = algo.evolve(pop)

    print('The optimum fitness value is {:0.12f} while the known/reference solution is {:0.12f}.'.format(pop.champion_f[0],refSol))               
    print('The optimum fitness occurs at {:0.12f} while the known/reference value is {:0.12f}.'.format(pop.champion_x[0],refX))               



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

