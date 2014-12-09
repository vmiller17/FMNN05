import numpy as N
import pylab as P 

from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem


#Bouncing ball
def f(t,y): 
    yd_1 = y[1]
    yd_2 = - 9.81
    return N.array([yd_1,yd_2])
    
    
    
    