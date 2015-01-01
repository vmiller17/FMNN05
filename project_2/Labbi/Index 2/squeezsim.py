# -*- coding: utf-8 -*-
from  __future__  import division
from  scipy       import *
import pylab as plt
import numpy as np

from assimulo.problem import Implicit_Problem
from assimulo.solvers import IDA
import squeezer

y0,yp0=squeezer.init_squeezer() # Initial values
t0 = 0 # Initial time

squeezemodel = Implicit_Problem(squeezer.squeezer, y0, yp0, t0)

algvar = 7*[1.]+13*[0.]
sim = IDA(squeezemodel) # Create solver instance
sim.algvar = algvar
sim.suppress_alg=True
sim.atol = 1e-7
tf = .03 # End time for simulation

t, y, yd = sim.simulate(tf)

fig = figure()
sim.plot(mask=7*[1]+13*[0])
plt.grid(1)
plt.axis([0, .03, -1.5, 1.5])
plt.xlabel('Time, t [s]')
plt.ylabel('Angle, [rad]')
plt.show()
fig.savefig('p2_index3_IDA.eps', format='eps',dpi=1000)


