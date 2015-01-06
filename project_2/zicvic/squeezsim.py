# -*- coding: utf-8 -*-
from  __future__  import division
from  scipy       import *
from  matplotlib.pyplot import *
import numpy as np

from assimulo.problem import Implicit_Problem
from assimulo.solvers import IDA
import squeezer_ind2
import squeezer_ind3

y20,yp20=squeezer_ind2.init_squeezer2() # Initial values
y30,yp30=squeezer_ind3.init_squeezer3() # Initial values
t0 = 0 # Initial time

squeezemodel2 = Implicit_Problem(squeezer_ind2.squeezer2, y20, yp20, t0)
squeezemodel3 = Implicit_Problem(squeezer_ind3.squeezer3, y30, yp30, t0)

algvar = 7*[1.]+13*[0.]
solver2 = IDA(squeezemodel2)    # Create solver instance
solver3 = IDA(squeezemodel3)    # Create solver instance
solver2.algvar = algvar
solver3.algvar = algvar

solver2.suppress_alg=True
solver3.suppress_alg=True
solver2.atol = 1e-7
solver3.atol = 1e-7
tf = .03 # End time for simulation

t2, y2, yd2 = solver2.simulate(tf)
t3, y3, yd3 = solver3.simulate(tf)

figure(1)
plot(t2,y2[:,14:])
figure(2)
plot(t3,y3[:,14:])
show()