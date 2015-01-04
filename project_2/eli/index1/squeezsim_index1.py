# -*- coding: utf-8 -*-
"""
Created on Sun Dec 14 12:17:52 2014
@author: Malin
"""
from  __future__  import division
from  scipy       import *
from  matplotlib.pyplot import *
import numpy as np

from assimulo.solvers import Dopri5
from assimulo.solvers import CVode
from assimulo.solvers import RungeKutta4
from assimulo.problem import Explicit_Problem
from assimulo.solvers import ExplicitEuler
import squeezer_index1

y0=squeezer_index1.init_squeezer() # Initial values
t0 = 0 # Initial time

squeezemodel = Explicit_Problem(squeezer_index1.squeezer_index1, y0,t0)

sim = RungeKutta4(squeezemodel) # Create solver instance
#sim.h = 0.01/100
tf = .03 # End time for simulation
sim.atol = 1e-7
sim.h = 1e-4
t,y=sim.simulate(tf)
sim.plot(mask=7*[1]+7*[0])
#plot(t,y[:,0:2])
#plot(t, y[:,2:7])
grid(1)
axis([0, .03, -1.5, 1.5])
xlabel('Time, t [s]')
ylabel('Angle, [rad]')
savefig('RungeKutta4.eps',format='eps',dpi = 1000)

