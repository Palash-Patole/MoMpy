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


###########################################
## Definition: PyGMO problem ##############   
###########################################

class pagmo_TwoburnTransfer:
    
    def fitness(self,x):
        print('nothing yet')
        
    def get_bounds(self):
        print('nothing yet')
   
#####################################
########### User inputs #############    
#####################################

# Parameters definition

visualize = 2 # 0 - no visualizations
              # 1 - [UPDATE THIS]
              
######################################
######## Computations ################
######################################


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

