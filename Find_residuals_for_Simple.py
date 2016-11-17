#!/usr/bin/python
"""
Created on Thu Nov 17 10:24:52 2016

@author: charles

The purpose of this program is to fit a DE to the time series data in order
to find potential gene regulators.
"""
import numpy
import sys
import Smooth_Data
import matplotlib.pyplot as plt

def parabola_derivative(x):
    return 2*x
def exponential(y):
    return y

def simple_model(X,y,y_basal = 0,y_decay = 0,X_weights = 0):
    """
    Simple model where rate of change for target gene y is a function of
    the current value of y and a linear combination or regulator genes x
    """
    return 0

def second_order_Runge_Kutta(X,y_0,T,der):
    """
    Also known as the midpoint method
    2nd order runge kutta aproximation of gene expression function for step sizes to be
    determine by 
    the time array, T
    the values of the regulator gene, X 2d array for values of regulators at times
    the initial value of the target gene, y_0
    a function representing the change in the target gene relative to time as a 
        function of the current value of the target gene and the current value of
        the regulator gene, as well as parameters
    """
    Y_Estimate = []
    Y_Estimate.append(y_0)
    for index,new_time in enumerate(T[1::]):
        step_size = new_time - T[index]
        y_midpoint = Y_Estimate[index] + (step_size/float(2))*der(Y_Estimate[index])
        #t_midpoint = float(X[index]) + (step_size/float(2))
        #print(t_midpoint)
        new_y = float(Y_Estimate[index]) + float(step_size*der(y_midpoint))
        Y_Estimate.append(new_y)
    return Y_Estimate
        



#fit_attempts = 100
#perturbation_data = Smooth_Data.parse_TSV(sys.argv[1])
if __name__ == "__main__":
    Y_Estimate = second_order_Runge_Kutta([0,1,2,3,4,5,6,7,8,9,10],1,[0,1,2,3,4,5,6,7,8,9,10],exponential)
    xp = numpy.linspace(0,10, 100)
    real = numpy.exp(xp)
    plt.plot([0,1,2,3,4,5,6,7,8,9,10],Y_Estimate,xp,real)
    plt.show()
