# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Laboratory

Date of creation: Tue Jun 30 10:24:25 2020

Version: Base 

Description:
        1. Optimization of the two-burn transfer problem from AstrodynamicsProblems/TransferOrbits
        2. Validation against the optimum solution found using grid search in AstrodynamicsProblems_ver0.py, problem = 1 [Local testing] - 6 algos
        3. Plot of best fitness vs # generations for all six alogs
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
import matplotlib.pyplot as plt


###########################################
## Definition: PyGMO problem ##############   
###########################################

class pygmo_TwoburnTransfer:
    
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
seedSet = [123, 7000, 22, 44050, 800]

# Alorithms of interest
Algos = {1:'Differential Evolution',2:'Artifical Bee Colony',3:'Improved Harmony Search', 4: 'Particle Swarm Optimization'}
Algos[5] = 'Simple Genetic Algorithm'
Algos[6] = "Corana’s Simulated Annealing"
visualize = 1 # 0 - no visualizations
              # 1 - Showing errors in obtaining optimum fitness with all algorithms
######################################
######## Computations ################
######################################
# Creating a storage variable
Results = np.zeros((len(Algos),6))              

# Setting the problem
problem = pg.problem(pygmo_TwoburnTransfer(para))
             
for i in range(1,len(Algos)+1):
    for seed in seedSet:
        if len(seedSet)==1:
            print('\n\n\n')
            print("******************************************************************")
            print('The selected algorithm is: ', Algos[i])
            print("******************************************************************")
        
        if i==1:
            algo = pg.algorithm(pg.de(nGen,seed=seed))
        elif i==2:
            algo = pg.algorithm(pg.bee_colony(nGen, limit = 20,seed=seed))
        elif i==3:
            algo = pg.algorithm(pg.ihs(gen=nGen,seed=seed))
        elif i==4:
            algo = pg.algorithm(pg.pso(gen=nGen,seed=seed))
        elif i==5:
            algo = pg.algorithm(pg.sga(gen=nGen,seed=seed))
        elif i==6: 
            algo = pg.algorithm(pg.simulated_annealing(seed=seed))
        else:
            pass
        
        pop = pg.population(problem,10)
        pop = algo.evolve(pop)
        
        Results[i-1,0] = pop.champion_f[0] # Best fitness of the population
        Results[i-1,1] = refSol # Reference solution
        Results[i-1,2] = 100 * (Results[i-1,0]-Results[i-1,1])/Results[i-1,1] # % Error
        Results[i-1,3] = pop.champion_x[0] # Best fitness occurs at
        Results[i-1,4] = refX # Reference solution
        Results[i-1,5] = 100 * (Results[i-1,3]-Results[i-1,4])/Results[i-1,4] # % Error
        
        if len(seedSet)==1:
             print('The optimum fitness value is {:0.12f} while the known/reference solution is {:0.12f}.'.format(pop.champion_f[0],refSol))               
             print('The optimum fitness occurs at {:0.12f} while the known/reference value is {:0.12f}.'.format(pop.champion_x[0],refX))               



#######################################
######### Visualization ###############
#######################################

if visualize == 0:
    print('No visualizations are requested by the user.')
elif visualize == 1:
    plt.figure()
    for count in range(1,len(Algos)+1):
        plt.semilogy(count,abs(Results[count-1,2]),'*',label=Algos[count])
        plt.ylabel('% error in the computed fitness value')
        title = "Performance of various optimization algorithms \nfor solving two-burn orbit transfer problem"
        plt.title(title)
        plt.legend()
        plt.grid(True)
#elif visualize ==2:
#    # [EDIT THIS]  
else:
    print('Unrecognized input for the visualization option.')

