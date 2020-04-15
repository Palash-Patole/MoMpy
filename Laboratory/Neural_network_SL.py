# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: Laboratory
Date of creation: Thu Apr 16 00:51:57 2020

Version: Base

Description:
    1. Script taken from SoloLearn to model a simple neuron, train and use it
    2. Actual output: (a+b)*2, when a and b are inputs
    3. Trains the network through number of iterations specified by user
"""

from numpy import array, random, dot


class neural_network:
    def __init__(self):
        random.seed(1)
        # we model a single neuron, with 3 inputs and 1 output and assign random weight        
        self.weights = 2 * random.random((2,1)) -1 # 2 by 1 array
        
    def train(self,inputs,outputs,num): # training the network for num times
        for iteration in range(num):
            output = self.think(inputs)
            error = outputs-output
            adjustment = 0.01 * dot(inputs.T, error) # adjustment = 0.01 * error * input
            self.weights += adjustment
            
    def think(self,inputs):
        return(dot(inputs,self.weights)) # output = weight1 * input1 + weight2 * input2
        
neural_network = neural_network()        

# The training set
inputs = array([[2,3], [1,1], [5,2], [12,3]])
outputs = array([[10, 4, 14, 30]]).T

# Training the neural network using the training set.
neural_network.train(inputs,outputs, 10000)


# Ask the neural network the output
print(neural_network.think(array([15,2])))        